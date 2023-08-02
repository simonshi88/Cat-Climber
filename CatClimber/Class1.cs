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



namespace CatClimber
{
    public class Script_Instance
    {
        private void RunScript(int sizeX, int sizeY, int sizeZ, bool reset, ref object A, ref object B)
        {
            if (!reset)
            {
                bool addConnections = true;

                PointGridManager gridManager = new PointGridManager(sizeX, sizeY, sizeZ);

                if (addConnections)
                {

                    gridManager.AddRandomConnections();
                    // 这里可以添加适当的逻辑来生成 0 和 1 的随机点阵，并调用 gridManager.AddConnection() 来添加连接。
                    // 假设我们这里简单地随机生成一个点阵，然后添加随机连接。
                    //        Random random = new Random();
                    //        for (int x = 0; x < sizeX; x++)
                    //        {
                    //          for (int y = 0; y < sizeY; y++)
                    //          {
                    //            for (int z = 0; z < sizeZ; z++)
                    //            {
                    //              int connectionChance = random.Next(10);
                    //              if (connectionChance < 5)
                    //              {
                    //                gridManager.AddConnection(x, y, z, random.Next(sizeX), random.Next(sizeY), random.Next(sizeZ));
                    //              }
                    //            }
                    //          }
                    //        }
                }


                points = new List<Point3d>();
                lines = new List<Line>();


                // 生成点和线
                for (int x = 0; x < sizeX; x++)
                {
                    for (int y = 0; y < sizeY; y++)
                    {
                        for (int z = 0; z < sizeZ; z++)
                        {
                            points.Add(new Point3d(x, y, z));

                            int connectionsX = gridManager.CalculateHorizontalConnections(x, y, z);
                            int connectionsY = gridManager.CalculateHorizontalConnectionsFB(x, y, z);
                            int connectionsZ = gridManager.CalculateVerticalConnections(x, y, z);
                            // 根据点与点之间的连接关系生成线
                            if (connectionsX > 0)
                            {
                                lines.Add(new Line(new Point3d(x, y, z), new Point3d(x + 1, y, z)));
                            }
                            if (connectionsY > 0)
                            {
                                lines.Add(new Line(new Point3d(x, y, z), new Point3d(x, y+1, z)));
                            }
                            if (connectionsZ > 0)
                            {
                                lines.Add(new Line(new Point3d(x, y, z), new Point3d(x, y, z + 1)));
                            }
                            // 在这里添加其他方向的线连接规则，根据实际情况调整。
                        }
                    }
                }

            }
            A = points;
            B = lines;

        }

        // <Custom additional code> 

        List<Point3d> points = new List<Point3d>();
        List<Line> lines = new List<Line>();
        public class PointGridManager
        {
            private int[,,] adjacencyMatrix; // 邻接矩阵
            private int gridSizeX; // 平面点阵的水平大小
            private int gridSizeY; // 平面点阵的垂直大小
            private int gridSizeZ; // 平面点阵的 z 方向大小

            private Random random;

            private List<Line> lines;

            public PointGridManager(int sizeX, int sizeY, int sizeZ)
            {
                gridSizeX = sizeX;
                gridSizeY = sizeY;
                gridSizeZ = sizeZ;
                adjacencyMatrix = new int[sizeX, sizeY, sizeZ];

                random = new Random();

                lines = new List<Line>();
            }

            public void ConnectPointsBasedOnAdjacencyMatrix()
            {
                for (int x = 0; x < gridSizeX; x++)
                {
                    for (int y = 0; y < gridSizeY; y++)
                    {
                        for (int z = 0; z < gridSizeZ; z++)
                        {
                            if (adjacencyMatrix[x, y, z] == 1)
                            {
                                List<Point3d> validNeighbors = GetValidNeighbors(x, y, z);
                                foreach (Point3d neighbor in validNeighbors)
                                {
                                    if (adjacencyMatrix[(int)Math.Round(neighbor.X), (int)Math.Round(neighbor.Y), (int)Math.Round(neighbor.Z)] == 1)
                                    {

                                        AddConnection(x, y, z, neighbor.X, neighbor.Y, neighbor.Z);
                                    }
                                }
                            }
                        }
                    }
                }
            }



            // 添加随机连接关系
            public void AddRandomConnections()
            {
                for (int x = 0; x < gridSizeX; x++)
                {
                    for (int y = 0; y < gridSizeY; y++)
                    {
                        for (int z = 0; z < gridSizeZ; z++)
                        {
                            int connectionChance = random.Next(10);
                            if (connectionChance < 5)
                            {
                                List<Point3d> validNeighbors = GetValidNeighbors(x, y, z);
                                if (validNeighbors.Count > 0)
                                {
                                    int randomNeighborIndex = random.Next(validNeighbors.Count);
                                    Point3d neighbor = validNeighbors[randomNeighborIndex];
                                    AddConnection(x, y, z, neighbor.X, neighbor.Y, neighbor.Z);
                                }
                            }
                        }
                    }
                }
            }

            // 获取可选的邻居点列表
            private List<Point3d> GetValidNeighbors(int x, int y, int z)
            {
                List<Point3d> validNeighbors = new List<Point3d>();

                if (IsValidNeighbor(x + 1, y, z))
                    validNeighbors.Add(new Point3d(x + 1, y, z));
                if (IsValidNeighbor(x - 1, y, z))
                    validNeighbors.Add(new Point3d(x - 1, y, z));
                if (IsValidNeighbor(x, y + 1, z))
                    validNeighbors.Add(new Point3d(x, y + 1, z));
                if (IsValidNeighbor(x, y - 1, z))
                    validNeighbors.Add(new Point3d(x, y - 1, z));
                if (IsValidNeighbor(x, y, z + 1))
                    validNeighbors.Add(new Point3d(x, y, z + 1));
                if (IsValidNeighbor(x, y, z - 1))
                    validNeighbors.Add(new Point3d(x, y, z - 1));

                return validNeighbors;
            }

            // 判断邻居点是否是有效的连接点
            private bool IsValidNeighbor(int x, int y, int z)
            {
                if (x < 0 || x >= gridSizeX || y < 0 || y >= gridSizeY || z < 0 || z >= gridSizeZ)
                    return false;

                // 在这里添加其他连接规则的判断逻辑
                // 例如：满足某些条件才能连线

                // 这里我们假设所有邻居都是有效的连接点
                return true;
            }

            // 添加连接关系
            public void AddConnection(int point1X, int point1Y, int point1Z, double neighborX, double neighborY, double neighborZ)
            {
                int point2X = (int)Math.Round(neighborX);
                int point2Y = (int)Math.Round(neighborY);
                int point2Z = (int)Math.Round(neighborZ);

                adjacencyMatrix[point1X, point1Y, point1Z] = 1;
                adjacencyMatrix[point2X, point2Y, point2Z] = 1;
                

            }

            // 计算水平方向连线数量
            public int CalculateHorizontalConnections(int pointX, int pointY, int pointZ)
            {
                int connectionsCount = 0;

                // 判断左右方向是否有连接
                if (pointX < gridSizeX - 1 && adjacencyMatrix[pointX + 1, pointY, pointZ] == 1)
                {
                    connectionsCount++;
                }
                if (pointX > 0 && adjacencyMatrix[pointX - 1, pointY, pointZ] == 1)
                {
                    connectionsCount++;
                }

                return connectionsCount;
            }

            public int CalculateHorizontalConnectionsFB(int pointX, int pointY, int pointZ)
            {
                int connectionsCount = 0;
                // 判断前后方向是否有连接
                if (pointY < gridSizeY - 1 && adjacencyMatrix[pointX, pointY + 1, pointZ] == 1)
                {
                    connectionsCount++;
                }
                if (pointY > 0 && adjacencyMatrix[pointX, pointY - 1, pointZ] == 1)
                {
                    connectionsCount++;
                }

                return connectionsCount;
            }

            // 计算垂直方向连线数量
            public int CalculateVerticalConnections(int pointX, int pointY, int pointZ)
            {
                int connectionsCount = 0;

                // 判断上下方向是否有连接
                if (pointZ < gridSizeZ - 1 && adjacencyMatrix[pointX, pointY, pointZ + 1] == 1)
                {
                    connectionsCount++;
                }
                if (pointZ > 0 && adjacencyMatrix[pointX, pointY, pointZ - 1] == 1)
                {
                    connectionsCount++;
                }

                return connectionsCount;
            }
        }


        // </Custom additional code> 
    }
}
