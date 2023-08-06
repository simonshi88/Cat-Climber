using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.DocObjects;


namespace CatClimber
{
    public class Script_Instance
    {
        private void RunScript(int sizeX, int sizeY, int sizeZ, bool reset, ref object A, ref object B)
        {
            if (!reset)
            {
                 manager = new SpaceConnectivityManager(sizeX, sizeY, sizeZ);

                // Randomly initialize connection matrices
                RandomlyInitializeMatrix(manager.ConnectionMatrixX);
                RandomlyInitializeMatrix(manager.ConnectionMatrixY);
                RandomlyInitializeMatrix(manager.ConnectionMatrixZ);

                // Get random lines based on the initialized matrices
                List<Line> randomLinesX = manager.GetConnectedLinesX();
                List<Line> randomLinesY = manager.GetConnectedLinesY();
                List<Line> randomLinesZ = manager.GetConnectedLinesZ();


                manager.GenerateLinesFromConnectionMatrix();

                All.AddRange(randomLinesX);
                All.AddRange(randomLinesY);
                All.AddRange(randomLinesZ);
                Points = manager.points;
            }
                

            A = manager.lines;
            B = Points ;
        }
        SpaceConnectivityManager manager;

        List<Point3d> Points = new List<Point3d>();
        List<Line> All = new List<Line>();

        // Helper function to randomly initialize a boolean matrix
        private static void RandomlyInitializeMatrix(bool[,,] matrix)
        {
            Random random = new Random();

            for (int x = 0; x < matrix.GetLength(0); x++)
            {
                for (int y = 0; y < matrix.GetLength(1); y++)
                {
                    for (int z = 0; z < matrix.GetLength(2); z++)
                    {
                        matrix[x, y, z] = random.Next(2) == 1; // 0 or 1
                    }
                }
            }
        }


        public class SpaceConnectivityManager
        {
            public bool[,,] ConnectionMatrixX;
            public bool[,,] ConnectionMatrixY;
            public bool[,,] ConnectionMatrixZ;
            public List<Point3d> points;
            public List<Line> lines;
            public SpaceConnectivityManager(int sizeX, int sizeY, int sizeZ)
            {
                ConnectionMatrixX = new bool[sizeX - 1, sizeY, sizeZ];
                ConnectionMatrixY = new bool[sizeX, sizeY - 1, sizeZ];
                ConnectionMatrixZ = new bool[sizeX, sizeY, sizeZ - 1];
                InitializePoints(sizeX, sizeY, sizeZ);
               
            }

            private void InitializePoints(int sizeX, int sizeY, int sizeZ)
            {
                points = new List<Point3d>();

                for (int x = 0; x < sizeX; x++)
                {
                    for (int y = 0; y < sizeY; y++)
                    {
                        for (int z = 0; z < sizeZ; z++)
                        {
                            points.Add(new Point3d(x, y, z));
                            
                        }
                    }
                }
            }

            public void GenerateLinesFromConnectionMatrix()
            {
                lines = new List<Line>();

                lines.AddRange(GenerateLinesFromConnectionMatrixX());
                lines.AddRange(GenerateLinesFromConnectionMatrixY());
                lines.AddRange(GenerateLinesFromConnectionMatrixZ());
               
            }

            public List<Line> GenerateLinesFromConnectionMatrixX()
            {
                List<Line> lines = new List<Line>();

                for (int x = 0; x < ConnectionMatrixX.GetLength(0); x++)
                {
                    for (int y = 0; y < ConnectionMatrixX.GetLength(1); y++)
                    {
                        for (int z = 0; z < ConnectionMatrixX.GetLength(2); z++)
                        {
                            if (ConnectionMatrixX[x, y, z])
                            {
                                Point3d startPoint = points[GetIndex(ConnectionMatrixX, x, y, z)];
                                Point3d endPoint = points[GetIndex(ConnectionMatrixX,x + 1, y, z)];
                                lines.Add(new Line(startPoint, endPoint));
                            }
                        }
                    }
                }

                return lines;
            }
            public List<Line> GenerateLinesFromConnectionMatrixY()
            {
                List<Line> lines = new List<Line>();

                for (int x = 0; x < ConnectionMatrixY.GetLength(0); x++)
                {
                    for (int y = 0; y < ConnectionMatrixY.GetLength(1); y++)
                    {
                        for (int z = 0; z < ConnectionMatrixY.GetLength(2); z++)
                        {
     
                            if (ConnectionMatrixY[x, y, z])
                            {
                                Point3d startPoint = points[GetIndex(ConnectionMatrixY, x, y, z)];
                                Point3d endPoint = points[GetIndex(ConnectionMatrixY, x, y + 1, z)];
                                lines.Add(new Line(startPoint, endPoint));
                            }

                        }
                    }
                }

                return lines;
            }

            public List<Line> FindLinesContainingPoint(Point3d targetPoint, List<Line> lines)
            {
                List<Line> containingLines = new List<Line>();

                foreach (Line line in lines)
                {

                    Point3d closestPoint = line.ClosestPoint(targetPoint,true);
                    if (closestPoint.DistanceTo(targetPoint) < Rhino.RhinoMath.SqrtEpsilon)
                    {
                        containingLines.Add(line);
                    }
                }

                return containingLines;
            }

            public DataTree<Line> FindLinesContainingPoint(DataTree<Point3d> points, List<Line> lines)
            {
                DataTree<Line> containingLines = new DataTree<Line>();

                for (int i = 0; i < points.BranchCount; i++)
                {
                    GH_Path path = points.Paths[i];
                    List<Point3d> pointBranch = points.Branch(path);
                    List<Line> lineBranch = new List<Line>();

                    foreach (Point3d point in pointBranch)
                    {
                        
                        foreach (Line line in lines)
                        {
                            Point3d closestPoint = line.ClosestPoint(point, true);

                            if (closestPoint.DistanceTo(point) < Rhino.RhinoMath.SqrtEpsilon)
                            {
                                lineBranch.Add(line);
                                   }

                        }
                        
                    }
                    containingLines.AddRange(lineBranch, path);

                }

                return containingLines;
            }

            public DataTree<Point3d> FindLayerPoints(List<Point3d> points)
            {
                DataTree<Point3d> layerPoints = new DataTree<Point3d>();

                for (int i = 0; i < points.Count; i++)
                {
                    int layer = (int)Math.Floor(points[i].Z);
                    layerPoints.Add(points[i], new GH_Path(layer));
                }
                return layerPoints;
            }


            public DataTree<Line> FindLayerLines(List<Line> line)
            {
                DataTree<Line> layerLines = new DataTree<Line>();

                for (int i = 0; i < line.Count; i++)
                {
                    if (line[i].From.Z == line[i].To.Z)
                    {
                        int layer = (int)Math.Floor(line[i].From.Z);
                        layerLines.Add(line[i], new GH_Path(layer));
            }

                }
                return layerLines;
            }


            public DataTree<Surface> FindLayerSurfaces(List<Surface> surfaces)
            {
                DataTree<Surface> layerSurfaces = new DataTree<Surface>();

                for (int i = 0; i < surfaces.Count; i++)
                {

                    int layer = (int)Math.Floor(surfaces[i].PointAt(0, 0).Z);
                    layerSurfaces.Add(surfaces[i], new GH_Path(layer));

                }
                return layerSurfaces;
            }


            public double CalculateMNGoal(List<int> m, List<int> n, double scoreStep)
            {
                double score = 0.0;

                if(m.Count < 1)
                    return score;

                if (m[0] > 2 && m[0] <= 4)
                {
                    if (n[0] == 4)
                    {
                        if (scoreStep > 0)
                            score += scoreStep;
                        else
                            score += 1;
                    }
                }
                else
                {
                    if (n[0] == m[0])
                    {
                        if (scoreStep > 0)
                            score += scoreStep;
                        else
                            score += 1;
                    }
                }

                for (int i = 1; i < m.Count; i++)
                {
                    if (m[i] >= 1 && m[i] < 4)
                    {
                        if (n[i] == m[i] + 1)
                        {
                            if (scoreStep > 0)
                                score += scoreStep;
                            else
                                score += 1;
                        }
                    }
                    else
                    {
                        if (n[i] == m[i])
                        {
                            if (scoreStep > 0)
                                score += scoreStep;
                            else
                                score += 1;
                        }
                    }
                }
                return score;
            }


            public DataTree<Point3d> CalculateNeighbors(List<Point3d> points, double spacing)
            {
                DataTree<Point3d> neighborsTree = new DataTree<Point3d>();

                for (int i = 0; i < points.Count; i++)
                {
                    List<Point3d> neighbors = new List<Point3d>();

                    for (int j = -1; j <= 1; j++)
                    {
                        for (int k = -1; k <= 1; k++)
                        {
                            for (int l = -1; l <= 1; l++)
                            {
                                if (j == 0 && k == 0 && l == 0)
                                    continue;

                                Point3d neighbor = new Point3d(
                                    points[i].X + j * spacing,
                                    points[i].Y + k * spacing,
                                    points[i].Z + l * spacing);

                                if (points.Contains(neighbor))
                                    neighbors.Add(neighbor);
                            }
                        }
                    }

                    neighborsTree.AddRange(neighbors, new GH_Path(i));
                }

                return neighborsTree;
            }


            public DataTree<int> CalculateConnectivityScoreAsDataTree(List<Point3d> points, List<Line> lines)
            {
                Dictionary<Point3d, int> connectivityScores = CalculateReachabilityCount(points, lines);

                DataTree<int> resultTree = new DataTree<int>();

                foreach (KeyValuePair<Point3d, int> kvp in connectivityScores)
                {
                    GH_Path path = new GH_Path(points.IndexOf(kvp.Key));
                    resultTree.Add(kvp.Value, path);
                }

                return resultTree;
            }


            public Dictionary<Point3d, int> CalculateReachabilityCount(List<Point3d> points, List<Line> lines)
            {
                Dictionary<Point3d, int> reachabilityCount = new Dictionary<Point3d, int>();

                // Create an adjacency list representation of the graph
                Dictionary<Point3d, List<Point3d>> adjacencyList = new Dictionary<Point3d, List<Point3d>>();
                foreach (Point3d point in points)
                {
                    adjacencyList[point] = new List<Point3d>();
                }

                foreach (Line line in lines)
                {
                    adjacencyList[line.From].Add(line.To);
                    adjacencyList[line.To].Add(line.From);
                }

                foreach (Point3d startPoint in points)
                {
                    HashSet<Point3d> visited = new HashSet<Point3d>();
                    Queue<Point3d> queue = new Queue<Point3d>();
                    queue.Enqueue(startPoint);
                    visited.Add(startPoint);

                    int reachableCount = 0;

                    while (queue.Count > 0)
                    {
                        Point3d currentPoint = queue.Dequeue();

                        foreach (Point3d neighbor in adjacencyList[currentPoint])
                        {
                            if (!visited.Contains(neighbor))
                            {
                                queue.Enqueue(neighbor);
                                visited.Add(neighbor);
                                reachableCount++;
                            }
                        }
                    }

                    reachabilityCount[startPoint] = reachableCount;
                }

                return reachabilityCount;
            }


            public DataTree<int> CalculateConnectivity( List<Line> lines, double distanceThreshold)
            {
                Dictionary<Point3d, int> connectivityScores = CalculateConnectivityScore(lines, distanceThreshold);

                DataTree<int> dataTree = new DataTree<int>();
                foreach (var kvp in connectivityScores)
                {
                    GH_Path path = new GH_Path((int)kvp.Key.X, (int)kvp.Key.Y, (int)kvp.Key.Z);
                    dataTree.Add(kvp.Value, path);
                }

                return dataTree;
            }


            public Dictionary<Point3d, int> CalculateConnectivityScore(List<Line> lines, double distanceThreshold)
            {
                // Create a dictionary to store line endpoints
                Dictionary<Point3d, List<Point3d>> lineEndpoints = new Dictionary<Point3d, List<Point3d>>();
                foreach (Line line in lines)
                {
                    if (!lineEndpoints.ContainsKey(line.From))
                        lineEndpoints[line.From] = new List<Point3d>();
                    lineEndpoints[line.From].Add(line.To);

                    if (!lineEndpoints.ContainsKey(line.To))
                        lineEndpoints[line.To] = new List<Point3d>();
                    lineEndpoints[line.To].Add(line.From);
                }

                // Create a dictionary to store connectivity scores
                Dictionary<Point3d, int> connectivityScores = new Dictionary<Point3d, int>();

                foreach (Point3d point in lineEndpoints.Keys)
                {
                    int connectivityScore = 0;

                    foreach (Point3d endpoint in lineEndpoints[point])
                    {
                        if (lineEndpoints.ContainsKey(endpoint) && lineEndpoints[endpoint].Contains(point))
                        {
                            connectivityScore++;
                        }
                    }

                    connectivityScores[point] = connectivityScore;
                }

                // Calculate local connectivity within 3x3x3 environment
                Dictionary<Point3d, int> localConnectivityScores = new Dictionary<Point3d, int>();

                foreach (Point3d point in lineEndpoints.Keys)
                {
                    int localConnectivityScore = 0;

                    foreach (Point3d nearbyPoint in lineEndpoints.Keys)
                    {
                        if (point.DistanceTo(nearbyPoint) <= distanceThreshold)
                        {
                            localConnectivityScore += connectivityScores.ContainsKey(nearbyPoint) ? connectivityScores[nearbyPoint] : 0;
                        }
                    }

                    localConnectivityScores[point] = localConnectivityScore;
                }

                return localConnectivityScores;
            }


            public void FindLines(List<Line> linesCollection,
                              out List<Line> allPlaneLines, out List<Line> allVerticalLines,
                              out DataTree<Line> xAxisParallelLines, out DataTree<Line> zAxisParallelLines, 
                              out DataTree<Line> layerPlaneLines,
                              out DataTree<Line> layerVerticalLines)
            {
                allPlaneLines = new List<Line>();
                allVerticalLines = new List<Line>();
                xAxisParallelLines = new DataTree<Line>();
                zAxisParallelLines = new DataTree<Line>();
                layerPlaneLines = new DataTree<Line>();
                layerVerticalLines = new DataTree<Line>();

                for (int i = 0; i < linesCollection.Count; i++)
                {
                    Line line = linesCollection[i];

                    if (line.Direction.IsParallelTo(Vector3d.ZAxis) != 0)
                    {
                        allVerticalLines.Add(line);
                    }
                    else
                    {
                        allPlaneLines.Add(line);
                    }


                    int layer = (int)Math.Floor(GetLineStartWithMinZ(line).Z);

                    layerPlaneLines.Add(line, new GH_Path(layer));
                    layerVerticalLines.Add(line, new GH_Path(layer));

                    if (line.Direction.IsParallelTo(Vector3d.XAxis) != 0)
                    {
                        xAxisParallelLines.Add(line, new GH_Path(layer));
                    }

                    if (line.Direction.IsParallelTo(Vector3d.ZAxis) != 0)
                    {
                        zAxisParallelLines.Add(line, new GH_Path(layer));
                    }
                }
            }



            public DataTree<double> ProcessZAxisParallelLines(DataTree<Line> zAxisParallelLinesByLayer, double arrayGap, int numberX, double scoreStep)
            {
                DataTree<double> allLayerGoal = new DataTree<double>();
                for (int i = 0; i < zAxisParallelLinesByLayer.BranchCount; i++)
                {
                    List<Line> lines = zAxisParallelLinesByLayer.Branch(i);
                    List<Point3d> startPoints = new List<Point3d>();

                    foreach (Line line in lines)
                    {
                        startPoints.Add(line.From);
                    }

                    Dictionary<double, List<Point3d>> groupedStartPoints = GroupStartPointsByY(startPoints);

                    List<double> goal = new List<double>();
                    foreach (var kvp in groupedStartPoints)
                    {                      
                        double roundedY = kvp.Key;
                        List<Point3d> points = kvp.Value;

                        double score = 0;
                        foreach (Point3d point in points)
                        {
                            for (int x = 0; x <= numberX - 2; x++)
                            {
                                List<Point3d> group = new List<Point3d> {
                                new Point3d(x * arrayGap, roundedY, points[0].Z),
                                    new Point3d((x+1) * arrayGap, roundedY, points[0].Z),
                                new Point3d((x+2) * arrayGap, roundedY, points[0].Z) };

                                
                                if (group.Contains(point))
                                {
                                    if (scoreStep > 0)
                                        score += scoreStep;
                                    else
                                        score += 1;                        
                                    
                                }

                            }
                            
                        }
                        goal.Add(score);
                    }

                    allLayerGoal.AddRange(goal, new GH_Path(i));
                }
                return allLayerGoal;

            }

            public Dictionary<double, List<Point3d>> GroupStartPointsByY(List<Point3d> startPoints)
            {
                Dictionary<double, List<Point3d>> groupedStartPoints = new Dictionary<double, List<Point3d>>();

                foreach (Point3d startPoint in startPoints)
                {
                    double roundedY = Math.Round(startPoint.Y, 2);

                    if (!groupedStartPoints.ContainsKey(roundedY))
                    {
                        groupedStartPoints[roundedY] = new List<Point3d>();
                    }

                    groupedStartPoints[roundedY].Add(startPoint);
                }

                return groupedStartPoints;
            }


            public double CalculateScore(double x, double mean, double stdDev, double minY, double maxY, double peakX)
            {
                // 计算正态分布曲线上的 y 值
                double y = NormalDistribution(x, mean, stdDev);

                // 根据 y 值的区间范围计算得分
                double score = MapToScore(y, minY, maxY);

                // 如果 x 等于指定的 peakX，则得分最大
                if (x == peakX)
                {
                    score = maxY;
                }

                return score;
            }

            // 计算正态分布函数
            private double NormalDistribution(double x, double mean, double stdDev)
            {
                double exponent = -(Math.Pow(x - mean, 2) / (2 * Math.Pow(stdDev, 2)));
                return Math.Exp(exponent) / (stdDev * Math.Sqrt(2 * Math.PI));
            }

            // 将 y 值映射到指定区间范围的得分
            private double MapToScore(double y, double minY, double maxY)
            {
                return minY + (maxY - minY) * y;
            }


            public Point3d GetLineStartWithMinZ(Line line)
            {
                if (line.From.Z <= line.To.Z)
                {
                    return line.From;
                }
                else
                {
                    return line.To;
                }
            }



            public void RemoveLinesBasedOnRandom(List<Line> lines)
            {
                Random random = new Random();
                List<Line> linesToRemove = new List<Line>();

                foreach (Line line in lines)
                {
                    int randomNumber = random.Next(0, 11); // Generates random number between 0 and 10

                    if (randomNumber < 5)
                    {
                        linesToRemove.Add(line);
                    }
                }

                foreach (Line lineToRemove in linesToRemove)
                {
                    lines.Remove(lineToRemove);
                }
            }

            public List<Line> FindVerticalLines(List<Line> lines)
            {
                List<Line> verticalLines = new List<Line>();

                foreach (Line line in lines)
                {
                    if (line.Direction.IsParallelTo(Vector3d.ZAxis) != 0)
                    {
                        verticalLines.Add(line);
                    }
                }

                return verticalLines;
            }


            public List<Line> GenerateLinesFromConnectionMatrixZ()
            {
                List<Line> lines = new List<Line>();

                for (int x = 0; x < ConnectionMatrixZ.GetLength(0); x++)
                {
                    for (int y = 0; y < ConnectionMatrixZ.GetLength(1); y++)
                    {
                        for (int z = 0; z < ConnectionMatrixZ.GetLength(2); z++)
                        {
      
                            if (ConnectionMatrixZ[x, y, z])
                            {
                                Point3d startPoint = points[GetIndex(ConnectionMatrixZ, x, y, z)];
                                Point3d endPoint = points[GetIndex(ConnectionMatrixZ, x, y, z + 1)];
                                lines.Add(new Line(startPoint, endPoint));
                            }
                        }
                    }
                }

                return lines;
            }






            public void AddConnectionX(int x, int y, int z)
            {
                if (IsValidPoint(ConnectionMatrixX,x, y, z))
                {
                    ConnectionMatrixX[x, y, z] = true;
                }
            }

            public void AddConnectionY(int x, int y, int z)
            {
                if (IsValidPoint(ConnectionMatrixY,x, y, z))
                {
                    ConnectionMatrixY[x, y, z] = true;
                }
            }

            public void AddConnectionZ(int x, int y, int z)
            {
                if (IsValidPoint(ConnectionMatrixZ,x, y, z))
                {
                    ConnectionMatrixZ[x, y, z] = true;
                }
            }

            public List<Line> GetConnectedLinesX()
            {
                return GetConnectedLinesFromMatrix(ConnectionMatrixX, new Vector3d(1, 0, 0));
            }

            public List<Line> GetConnectedLinesY()
            {
                return GetConnectedLinesFromMatrix(ConnectionMatrixY, new Vector3d(0, 1, 0));
            }

            public List<Line> GetConnectedLinesZ()
            {
                return GetConnectedLinesFromMatrix(ConnectionMatrixZ, new Vector3d(0, 0, 1));
            }

            private List<Line> GetConnectedLinesFromMatrix(bool[,,] matrix, Vector3d direction)
            {
                List<Line> lines = new List<Line>();

                for (int x = 0; x < matrix.GetLength(0); x++)
                {
                    for (int y = 0; y < matrix.GetLength(1); y++)
                    {
                        for (int z = 0; z < matrix.GetLength(2); z++)
                        {
                            if (matrix[x, y, z])
                            {
                                Point3d startPoint = points[GetIndex(matrix,x, y, z)];
                                Point3d endPoint = startPoint + direction;
                                lines.Add(new Line(startPoint, endPoint));
                            }
                        }
                    }
                }

                return lines;
            }

            private int GetIndex(bool[,,] matrix, int x, int y, int z)
            {
                return x + y * (matrix.GetLength(0)) + z * (matrix.GetLength(0) * matrix.GetLength(1));
            }

            private bool IsValidPoint(bool[,,] matrix, int x, int y, int z)
            {
                return x >= 0 && x < matrix.GetLength(0) &&
                       y >= 0 && y < matrix.GetLength(1) &&
                       z >= 0 && z < matrix.GetLength(2);
            }
        }





    }




    public class NeighborInfo
    {
        public Point3d Point { get; set; }
        public bool[] HasNeighborInDirection  = new bool[6]; // 0: X+, 1: X-, 2: Y+, 3: Y-, 4: Z+, 5: Z-
        public bool HasNeighborUp { get; set; }
        public bool HasNeighborDown { get; set; }
    }
    
    public class NeighborChecker
    {
        public DataTree<int> CalculateScores(List<Line> lines)
        {
            DataTree<int> result = new DataTree<int>();
            List<NeighborInfo> infoList = CheckNeighborsAndVerticalLines(lines);

            foreach (NeighborInfo info in infoList)
            {
                GH_Path path = new GH_Path((int)info.Point.X, (int)info.Point.Y, (int)info.Point.Z);

                int score = 0;
                if (info.HasNeighborInDirection[0] || info.HasNeighborInDirection[1] ||
                    info.HasNeighborInDirection[2] || info.HasNeighborInDirection[3] ||
                    info.HasNeighborInDirection[4] || info.HasNeighborInDirection[5])
                {
                    score += 1;
                }
                if (info.HasNeighborDown)
                {
                    score += 2;
                }
                if (info.HasNeighborUp)
                {
                    score += 2;
                }

                result.Add(score, path);
            }

            return result;
        }


        public List<NeighborInfo> CheckNeighborsAndVerticalLines(List<Line> lines)
        {
            List<NeighborInfo> result = new List<NeighborInfo>();

            foreach (Line line in lines)
            {
                // Check each endpoint of the line
                foreach (Point3d endpoint in new Point3d[] { line.From, line.To })
                {
                    NeighborInfo info = new NeighborInfo();
                    info.Point = endpoint;

                    for (int i = 0; i < 6; i++)
                    {
                        info.HasNeighborInDirection[i] = HasNeighborInDirection(endpoint, lines, (Direction)i);
                    }

                    info.HasNeighborUp = HasNeighborUp(endpoint, lines);
                    info.HasNeighborDown = HasNeighborDown(endpoint, lines);

                    result.Add(info);
                }
            }

            return result;
        }

        private bool HasNeighborInDirection(Point3d point, List<Line> lines, Direction direction)
        {
            Vector3d vector = GetDirectionVector(direction);

            foreach (Line neighborLine in lines)
            {
                if (neighborLine.From.DistanceTo(point) < 1e-6 || neighborLine.To.DistanceTo(point) < 1e-6)
                {
                    Vector3d neighborDirection = neighborLine.Direction;
                    if (vector.IsParallelTo(neighborDirection, 1e-6) != 0)
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        private bool HasNeighborUp(Point3d point, List<Line> lines)
        {
            foreach (Line neighborLine in lines)
            {
                if (neighborLine.From.DistanceTo(point) < 1e-6 || neighborLine.To.DistanceTo(point) < 1e-6)
                {
                    if (neighborLine.Direction.Z > 1e-6)
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        private bool HasNeighborDown(Point3d point, List<Line> lines)
        {
            foreach (Line neighborLine in lines)
            {
                if (neighborLine.From.DistanceTo(point) < 1e-6 || neighborLine.To.DistanceTo(point) < 1e-6)
                {
                    if (neighborLine.Direction.Z < -1e-6)
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        private Vector3d GetDirectionVector(Direction direction)
        {
            switch (direction)
            {
                case Direction.XPositive:
                    return new Vector3d(1, 0, 0);
                case Direction.XNegative:
                    return new Vector3d(-1, 0, 0);
                case Direction.YPositive:
                    return new Vector3d(0, 1, 0);
                case Direction.YNegative:
                    return new Vector3d(0, -1, 0);
                case Direction.ZPositive:
                    return new Vector3d(0, 0, 1);
                case Direction.ZNegative:
                    return new Vector3d(0, 0, -1);
                default:
                    return Vector3d.Zero;
            }
        }

        private enum Direction
        {
            XPositive,
            XNegative,
            YPositive,
            YNegative,
            ZPositive,
            ZNegative
        }
    }
}
