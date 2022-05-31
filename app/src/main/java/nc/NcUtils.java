package nc;

import android.util.Log;

public class NcUtils {

    static{
        //加载so依赖库
        System.loadLibrary("native-lib");;
    }

    public NcUtils() {
    }

    //申请GF(2^m)
    public static native void initGalois(int m);

    //释放jni申请的空间
    public static native void uninitGalois();

    //随机生成一个有限域上的数
    public static native int gf();

    //有限域加
    public static native int Add(int a,int b);

    //有限域减
    public static native int Sub(int a,int b);

    //有限域乘
    public static native int Mul(int a,int b);

    //有限域除
    public static native int Div(int a,int b);

    //有限域指数
    public static native int Exp(int a,int n);

    //有限域求逆元 即-1次
    public static native int Inv(int a);

    //生成有限域上随机矩阵
    public static native int[][] Matrix(int row,int col);

    //矩阵相乘
    public static native int[][] MatrixMultiply(int[][] matrix1,int[][] matrix2);

    //矩阵求秩
    public static native int Rank(int[][] matrix);

    //矩阵求逆
    public static native int[][] Inverse(int[][] matrix);

    //LU矩阵解码
    public static native int[][] InverseLU(int[][] matrix1,int[][] matrix2, int threads);

    //生成有限域上的随机稀疏矩阵
    public static native int[][] SparseMatrix(int n,int row,int col);

    //求方阵的行列式的值Det
    public static native int Det(int[][] matrix);

}

