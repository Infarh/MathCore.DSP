namespace MathCore.DSP;

internal static class MatrixZ
{
    public static void Test()
    {
        var Z = new int[6, 6];
        var (n, m) = Z;

        Z[0, 0] = 1;

        //  |  0   1   2   3   4   5 -> j
        // -+-----------------------
        // 0|  1   1   1   1   1   1
        // 1|      1   2   3   4   5
        // 2|          1   3   6  10
        // 3|              1   4  10
        // 4|                  1   5
        // 5|                      1
        // \|/
        //  i

        for (var j = 1; j < m; j++)
        {
            Z[0, j] = 1;
            Z[j, j] = 1;

            for (var i = 1; i < j; i++)
            {
                Z[i, j] = Z[i, j - 1] + Z[i - 1, j - 1];
            }
        }

        Z.Print();


        //  |  0   1   2   3   4   5 -> j
        // -+-----------------------
        // 0|  1   1   1   1   1   1
        // 1|      1   2   3   4   5
        // 2|          1   3   6  10
        // 3|              1   4  10
        // 4|                  1   5
        // 5|                      1
        // \|/
        //  i

        //  1     =  1
        //  4 -1  =  3
        //  6 -4  =  2
        //  4 -6  = -2
        //  1 -4  = -3
        //    -1  = -1

        //  1         =  1
        //  3 -2*1    =  1
        //  3 -2*3 1  = -2
        //  1 -2*3 3  = -2
        //    -2*1 3  =  1
        //         1  =  1

        //  |  0   1   2   3   4   5 -> j
        // -+-----------------------
        // 0|  1   1   1   1   1   1
        // 1|      1   2   3   4   5
        // 2|          1   3   6  10
        // 3|              1   4  10
        // 4|                  1   5
        // 5|                      1
        // \|/
        //  i

        Console.WriteLine();
        //for (var j = m - 2; j >= m / 2; j--)
        //{
        //    for(var i0 = 1; )

        //    break;
        //}

        Z.Print();
    }

    private static void Print(this int[,] M)
    {
        var (n, m) = M;

        for (var i = 0; i < n; i++)
        {
            for (var j = 0; j < m; j++)
            {
                if (M[i, j] == 0)
                    Console.Write("   .");
                else
                    Console.Write($"{M[i, j],4}");
            }
            Console.WriteLine();

        }
    }
}
