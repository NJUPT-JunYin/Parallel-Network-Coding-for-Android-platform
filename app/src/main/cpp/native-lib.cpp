#include <jni.h>
#include <string>
#include "stdlib.h"
#include <malloc.h>
#include <time.h>
#include <android/log.h>

//#include<string.h>
//#include<time.h>
//#include<iostream>
//#include<stdio.h>
//#include<stdlib.h>
//#include<malloc.h>
#include <omp.h>

//有限域的加减乘指数运算宏定义
#define  gf_alpha(n)   (table_alpha[n])
#define  gf_index(n)   (table_index[n])


#define  gf_add(a,b)	(a^b)
#define  gf_sub(a,b)	(a^b)
#define  gf_xor(a,b)    (a^b)

#define  gf_mul(a,b)	(table_mul[a][b])
#define  gf_div(a,b)	(table_div[a][b])

//加法负元与乘法逆元
#define  gf_neg(a)		(a)
#define  gf_inv(a)		(table_div[1][a])

typedef unsigned int GFType;
typedef unsigned char byte;

//有限域
int gFieldSize;

//乘法表和除法表
GFType *table_alpha;
GFType *table_index;
GFType **table_mul;
GFType **table_div;

//有限域表（修改后与matlab一致）
GFType prim_poly[13] ={
/*	0 */	0x00000000,
/*  1 */    0x00000001,
/*  2 */    0x00000007,
/*  3 */    0x0000000b,
/*  4 */    0x00000013,
/*  5 */    0x00000025,
/*  6 */    0x00000043,
/*  7 */    0x00000089,
///*  8 */    0x00000187,
              0x0000011d,
/*  9 */    0x00000211,
/* 10 */    0x00000409,
/* 11 */    0x00000805,
/* 12 */    0x00001053,
};

extern "C" {

//初始化有限域
void gf_init(int m, GFType prim);
//释放有限域
void gf_uninit();

//按公式计算，L、U矩阵中i行j列的累加和（K个累加和）
//jint N;


GFType gfexp(GFType a, GFType n);

GFType gfmul(GFType a, GFType b);

GFType gfdiv(GFType a, GFType b);

GFType gfxor(GFType a,GFType b);

GFType gf();

//int N;
//按公式计算，L、U矩阵中i行j列的累加和（K个累加和）
JNIEXPORT jint JNICALL
sumLU(int i, int j, int K, jint **L, jint **U,int n);

//初始化有限域（二维查表法）
JNIEXPORT void JNICALL
Java_nc_NcUtils_initGalois(JNIEnv *env,jclass ,jint m);

//释放有限域
JNIEXPORT void JNICALL
Java_nc_NcUtils_uninitGalois(JNIEnv *env,jclass );

//生成有限域上随机数
JNIEXPORT jint JNICALL
Java_nc_NcUtils_gf(JNIEnv *env,jclass );

//有限域加
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Add(JNIEnv *env,jclass ,jint a,jint b);

//有限域减
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Sub(JNIEnv *env,jclass ,jint a,jint b);

//有限域乘
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Mul(JNIEnv *env,jclass ,jint a,jint b);

//有限域除
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Div(JNIEnv *env,jclass ,jint a,jint b);

//有限域指数
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Exp(JNIEnv *env,jclass ,jint a,jint n);

//求逆元
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Inv(JNIEnv *env,jclass ,jint a);

//生成有限域随机矩阵
JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_Matrix(JNIEnv *env,jclass ,jint row,jint col);

//矩阵相乘
JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_MatrixMultiply(JNIEnv *env,jclass ,jobjectArray matrix1,jobjectArray matrix2);

//矩阵求秩
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Rank(JNIEnv *env,jclass ,jobjectArray matrix);

//矩阵求逆
JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_Inverse(JNIEnv *env,jclass ,jobjectArray matrix);


//矩阵求逆
JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_InverseLU(JNIEnv *env,jclass ,jobjectArray matrix1,jobjectArray matrix2,jint threads);

//有限域求方阵的行列式的值
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Det(JNIEnv *env,jclass ,jobjectArray matrix);

//生成有限域上的随机稀疏矩阵
JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_SparseMatrix(JNIEnv *env,jclass ,jint n,jint row,jint col);


}

void gf_init(int m, GFType prim){
    int i=0,j=0;

    srand((unsigned)(time(NULL)));      //获取时间种子 为随机数做准备

    if (m > 12 || m < 0)	// the field size is supported from GF(2^1) to GF(2^12).
        return;

    gFieldSize = 1<<m;

    if (0 == prim)
        prim = prim_poly[m];


    table_alpha = (GFType*)malloc(sizeof(GFType)*gFieldSize);
    table_index = (GFType*)malloc(sizeof(GFType)*gFieldSize);
    table_mul = (GFType**)malloc(sizeof(GFType*)*gFieldSize);
    table_div = (GFType**)malloc(sizeof(GFType*)*gFieldSize);
    for(i=0; i<gFieldSize; i++)
    {
        table_mul[i] = (GFType *)malloc(sizeof(GFType) * gFieldSize);
        table_div[i] = (GFType *)malloc(sizeof(GFType) * gFieldSize);
    }

    table_alpha[0]=1;
    table_index[0]=-1;

    for (i=1; i<gFieldSize; i++)
    {
        table_alpha[i] = table_alpha[i-1]<<1;
        if (table_alpha[i]>=gFieldSize)
        {
            table_alpha[i]^=prim;
        }
        table_index[table_alpha[i]]=i;
    }

    table_index[1]=0;

    // 生成乘法表与除法表
    for (i=0; i<gFieldSize; i++)
        for (j=0; j<gFieldSize; j++)
        {
            table_mul[i][j]=gfmul(i,j);
            table_div[i][j]=gfdiv(i,j);
        }
}


void gf_uninit(){
    int i = 0;

    free(table_alpha);
    free(table_index);

    for(i=0; i<gFieldSize; i++)
    {
        free(table_mul[i]);
        free(table_div[i]);
    }
    free(table_mul);
    free(table_div);
}

GFType gfexp(GFType a, GFType n){
    return table_alpha[table_index[a]*n%(gFieldSize-1)];
}

GFType gfmul(GFType a, GFType b) {
    if (0 == a || 0 == b)
        return 0;

    return table_alpha[(table_index[a] + table_index[b]) % (gFieldSize - 1)];
}

GFType gfdiv(GFType a, GFType b) {
    if (0 == a || 0 == b)
        return 0;

    return table_alpha[(table_index[a] - table_index[b] + (gFieldSize - 1)) % (gFieldSize - 1)];
}

GFType gfxor(GFType a,GFType b) {
    return table_alpha[(table_index[a]^table_index[b])];
}

GFType gf(){
    GFType i;
    i = (GFType)rand() % gFieldSize;
    return i;
}


JNIEXPORT void JNICALL
Java_nc_NcUtils_initGalois(JNIEnv *env,jclass ,jint m){
    gf_init(m,prim_poly[m]);
}

JNIEXPORT void JNICALL
Java_nc_NcUtils_uninitGalois(JNIEnv *env,jclass ){
    gf_uninit();
}

JNIEXPORT jint JNICALL
Java_nc_NcUtils_gf(JNIEnv *env,jclass ){
    return gf();
}

JNIEXPORT jint JNICALL
Java_nc_NcUtils_Add(JNIEnv *env,jclass ,jint a,jint b){
    return gf_add(GFType(a),GFType(b));
}

JNIEXPORT jint JNICALL
Java_nc_NcUtils_Sub(JNIEnv *env,jclass ,jint a,jint b){
    return gf_sub(GFType(a),GFType(b));
}

JNIEXPORT jint JNICALL
Java_nc_NcUtils_Mul(JNIEnv *env,jclass ,jint a,jint b){
    return gfmul(GFType(a),GFType(b));
}

JNIEXPORT jint JNICALL
Java_nc_NcUtils_Div(JNIEnv *env,jclass ,jint a,jint b){
    return gfdiv(GFType(a),GFType(b));
}

JNIEXPORT jint JNICALL
Java_nc_NcUtils_Exp(JNIEnv *env,jclass ,jint a,jint n){
    return gfexp(GFType(a),GFType(n));
}

JNIEXPORT jint JNICALL
Java_nc_NcUtils_Inv(JNIEnv *env,jclass ,jint a){
    return gf_inv(a);
}

JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_Matrix(JNIEnv *env,jclass ,jint row,jint col){
    //生成一个有限域上row*col的矩阵
    jint ** mat;
    int i,j;
    mat = new jint *[row];
    for(i=0;i<row;i++){
        mat[i] = new jint[col];
    }
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            mat[i][j] = (jint) gf();
        }
    }

    //将矩阵传回java层
    jclass intClass =env->FindClass("[I");
    jobjectArray matrix = env->NewObjectArray(row,intClass,NULL);
    for(i=0;i<row;i++){
        jintArray intdata = env->NewIntArray(col);
        env->SetIntArrayRegion(intdata,0,col,*&mat[i]);
        env->SetObjectArrayElement(matrix,i,intdata);
        env->DeleteLocalRef(intdata);
    }
    return matrix;
}

JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_MatrixMultiply(JNIEnv* env,jclass ,jobjectArray matrix1,jobjectArray matrix2){
    int i,j,k;

    //获取矩阵1的行列数
    int row1 = (int)(env->GetArrayLength(matrix1));
    jarray arr1 = (jarray)(env->GetObjectArrayElement(matrix1,0));
    int col1 = (int)(env->GetArrayLength(arr1));

    //获取矩阵2的行列数
    int row2 = (int)(env->GetArrayLength(matrix2));
    jarray arr2 = (jarray)(env->GetObjectArrayElement(matrix2,0));
    int col2 = (int)(env->GetArrayLength(arr2));

    //判断是否满足可乘条件
    if(col1 != row2){
        return NULL;
    }

    //获取矩阵1
    jint mat1[row1][col1];
    for(i=0;i<row1;i++){
        arr1 = (jarray)(env->GetObjectArrayElement(matrix1,i));
        jint *coldata = (jint*)(env->GetIntArrayElements((jintArray)arr1,0));
        for(j=0;j<col1;j++){
//            for(k=0;k<row2;k++) //
            mat1[i][j] = coldata[j];
        }
        env->ReleaseIntArrayElements((jintArray)arr1,coldata,0);
    }

    //获取矩阵2
    jint mat2[row2][col2];
    for(i=0;i<row2;i++){
        arr2 = (jarray)(env->GetObjectArrayElement(matrix2,i));
        jint *coldata = (jint*)(env->GetIntArrayElements((jintArray)arr2,0));
        for(j=0;j<col2;j++){
//            for(k=0;k<row2;k++) //
            mat2[i][j] = coldata[j];
        }
        env->ReleaseIntArrayElements((jintArray)arr2,coldata,0);
    }

    //矩阵相乘相关代码
    GFType temp;
    jint ** mat3;
    mat3 = new jint *[row1];
    for(i=0;i<row1;i++){
//        for(k=0;k<row2;k++) //
        mat3[i] = new jint[col2];
    }
    for(i=0;i<row1;i++){
        for(j=0;j<col2;j++){
            temp=0;
            for(k=0;k<row2;k++){
//                for(int q=0;q<row2;q++) //
                temp = gf_add(temp,gf_mul(mat1[i][k],mat2[k][j]));
            }
            mat3[i][j] = (GFType) temp;
        }
    }

    //返回结果
    jclass intClass =env->FindClass("[I");
    jobjectArray matrix3 = env->NewObjectArray(row1,intClass,NULL);
    for(i=0;i<row1;i++){
        jintArray intdata = env->NewIntArray(col2);
        env->SetIntArrayRegion(intdata,0,col2,*&mat3[i]);
        env->SetObjectArrayElement(matrix3,i,intdata);
        env->DeleteLocalRef(intdata);
    }
    return matrix3;
}

JNIEXPORT jint JNICALL
Java_nc_NcUtils_Rank(JNIEnv *env,jclass ,jobjectArray matrix){
    int i,j;

    //获取矩阵行列
    int row = (int)(env->GetArrayLength(matrix));
    jarray arr = (jarray)(env->GetObjectArrayElement(matrix,0));
    int col = (int)(env->GetArrayLength(arr));

    //获取矩阵
    jfloat mat[row][col];
    for(i=0;i<row;i++){
        arr = (jarray)(env->GetObjectArrayElement(matrix,i));
        jint *coldata = (jint *)(env->GetIntArrayElements((jintArray)arr,0));
        for(j=0;j<col;j++){
            mat[i][j] = coldata[j];
        }
        env->ReleaseIntArrayElements((jintArray)arr,coldata,0);
    }

    int m=0,n=0,k,r,u;
    jfloat temp;
    jfloat div;
    while(m<row && n<col){
        r = m;      //选主行
        for(k=m;k<row;k++){
            if(mat[k][n]){
                r=k;
                break;
            }
        }
        if (mat[r][n]){
            if( r != m){
                for(k=0;k<=col;k++){
                    temp = mat[r][k];
                    mat[r][k] = mat[m][k];
                    mat[m][k] = temp;
                }
            }
            for(u=m+1;u<row;u++){
                if(mat[u][n]){
                    div = mat[u][n] / mat[m][n];
                    for(k=m;k<col;k++){
                        mat[u][k] -= div * mat[m][k];
                    }
                }
            }
            m++;
        }
        n++;
    }
    return m;
}
//高斯消元函数！-------------------------------------------------------------
//矩阵求逆
JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_Inverse(JNIEnv *env,jclass jclass1,jobjectArray matrix){
    int i,j,x,z,k;

    //获取矩阵行列
    int row = (int)(env->GetArrayLength(matrix)); //调用函数 得到矩阵长度----行
    jarray arr = (jarray)(env->GetObjectArrayElement(matrix,0));//调用函数 得到矩阵第一行元素构成的矩阵？
    int col = (int)(env->GetArrayLength(arr));    //调用函数 得到矩阵长度----列

    //不是方阵返回null
    if(row != col || row != (int)Java_nc_NcUtils_Rank(env,jclass1,matrix)){
        return NULL;
    }

    //获取增广矩阵
    jint mat[row][row*2]; //构造矩阵   列数变为2倍  为单位矩阵设置空间
    for(i=0;i<row;i++){   //前row列    为原始矩阵
        arr = (jarray)(env->GetObjectArrayElement(matrix,i)); //原始矩阵元素第i行存放到arr
        jint *coldata = (jint *)(env->GetIntArrayElements((jintArray)arr,0));//arr 转化到coldata? 为什么？
        for(j=0;j<row;j++){
            mat[i][j] = coldata[j];   //元素值依次赋值到新矩阵中
        }
        env->ReleaseIntArrayElements((jintArray)arr,coldata,0);  //释放？ 释放空间？
        for(j=row;j<row*2;j++){  //后row列    构造单位矩阵  对角线为1
            if(i == j - row){
                mat[i][j] = 1;
            }else{
                mat[i][j] = 0;
            }
        }
    }

    //step1:将增广矩阵转换成下三角矩阵
    for (i = 0; i < row; i++) {
        if (mat[i][i] == 0) {     //如果出现主对角线元素值为0的情况  与其他行进行交换
            for (z = i + 1; z < row; z++) {    //下一行相同列元素  不为0就开始交换
                if (mat[z][i] != 0) {
                    for (x = 0; x < 2 *row; x++) {  //两行的元素值  依次交换
                        int temp = mat[i][x];
                        mat[i][x] = mat[z][x];
                        mat[z][x] = temp;
                    }
                    break;
                }
            }
        }
        for (j = i + 1; j < row; j++)
        {   //按消元法思想  通过除法 得到系数
            GFType temp = gf_div(mat[j][i], mat[i][i]);
            for (z = 0; z < 2 * row; z++){   //根据系数  每一列进行初等变换  元素值变为0
                for (k = 2 * row; k >= 0; k--)
                    mat[j][z] = gf_add(mat[j][z], gf_mul(temp, mat[i][z]));
            }
        }

    }
    //step2:高斯消元
    for (i = 1; i < row; i++)
    {
        for (k = 0; k < i; k++)
        {  //按消元法思想  通过除法 得到系数
            GFType temp = gf_div(mat[k][i], mat[i][i]);
            for (z = 0; z < 2 * row; z++){   //根据系数  每一列进行初等变换  元素值变为0
                for (j = 2 * row; j >= 0; j--)
                    mat[k][z] = gf_add(mat[k][z], gf_mul(temp, mat[i][z]));
            }
        }
    }
    for (i = 0; i < row; i++) {
        if (mat[i][i] != 1)
        {    //主对角线变为1  构造单位矩阵
            jint temp = (GFType)mat[i][i];
            for (z = 0; z < 2 * row; z++){  //根据系数  每一列进行初等变换  元素值变为1
                for (j = 2 * row; j >= 0; j--){
                    mat[i][z] = gf_div(mat[i][z], temp);
                }
            }
        }
    }
//    int n_threads = 4;
//    omp_set_num_threads(n_threads);//设置线程数-----------------------------
    //获取求逆得到的矩阵
    jint result[row][row];   //构造结果矩阵 result
    for(i=0;i<row;i++)
    {
//#pragma omp parallel for   //并行化-----------------------------
        for(j=0;j<row;j++)
        {    //将第row+1~2row行元素输出到新的结果矩阵
            result[i][j] = mat[i][j+row];
        }
    }

    //返回结果                  //为什么还要另外构造一个新的矩阵？ JAVA用法决定？？？
    jclass intClass =env->FindClass("[I");
    jobjectArray mat_inv = env->NewObjectArray(row,intClass,NULL);
    for(i=0;i<row;i++){
        jintArray intdata = env->NewIntArray(row);
        env->SetIntArrayRegion(intdata,0,row,*&result[i]);
        env->SetObjectArrayElement(mat_inv,i,intdata);
        env->DeleteLocalRef(intdata);
    }
    return mat_inv;
}
//高斯消元函数！-------------------------------------------------------------

//LU运算-------------------------------------------------------------
//按公式计算，L、U矩阵中i行j列的累加和（K个累加和）
JNIEXPORT jint JNICALL
sumLU(int i, int j, int K, jint **L, jint **U,int n)
{
    GFType res = 0;
    for (int k = 0; k < K; k++)
    {
//        res += L[i][k] * U[k][j];
        res = gf_add(res,gf_mul(*((jint *)L+n*i+k), *((jint *)U+n*k+j)));//累加
    }
    return res;
}
//LU运算-------------------------------------------------------------


//////LU串行分解函数！--------------------------------------------------------------
//////矩阵求逆
//
//JNIEXPORT jobjectArray JNICALL
//Java_nc_NcUtils_InverseLU(JNIEnv *env,jclass jclass1,jobjectArray matrix1,jobjectArray matrix2){
//    int n;
//    int k,i,j;
//    int flag=1;
//    jint s,t;
//    //线程数----------------------------------------
////    int n_threads = 2;
//
//    //获取矩阵1的行列数
//    int row1 = (int)(env->GetArrayLength(matrix1));//调用函数 得到矩阵长度----行
//    jarray arr1 = (jarray)(env->GetObjectArrayElement(matrix1,0));//调用函数 得到矩阵第一行元素构成的矩阵？
//    int col1 = (int)(env->GetArrayLength(arr1));//调用函数 得到矩阵长度----列
//
//    //获取矩阵2的行列数
//    int row2 = (int)(env->GetArrayLength(matrix2));
//    jarray arr2 = (jarray)(env->GetObjectArrayElement(matrix2,0));
//    int col2 = (int)(env->GetArrayLength(arr2));
//
//    //判断是否满足可乘条件
//    if(col1 != row2){
//        return NULL;
//    }
//    //不是方阵返回null
//    if(row1 != col1 || row1 != (int)Java_nc_NcUtils_Rank(env,jclass1,matrix1)){
//        return NULL;
//    }
//    //将输入矩阵转化为函数中矩阵
//    n=row1;
//
//    //获取矩阵1
//    jint a[row1][col1];
//    for(i=0;i<row1;i++){
//        arr1 = (jarray)(env->GetObjectArrayElement(matrix1,i));
//        jint *coldata = (jint*)(env->GetIntArrayElements((jintArray)arr1,0));
//        for(j=0;j<col1;j++){
//            a[i][j] = coldata[j];
//        }
//        env->ReleaseIntArrayElements((jintArray)arr1,coldata,0);
//    }
//
//    //获取矩阵2
//    jint pDatMat[row2][col2];
//    for(i=0;i<row2;i++){
//        arr2 = (jarray)(env->GetObjectArrayElement(matrix2,i));
//        jint *coldata = (jint*)(env->GetIntArrayElements((jintArray)arr2,0));
//        for(j=0;j<col2;j++){
//            pDatMat[i][j] = coldata[j];
//        }
//        env->ReleaseIntArrayElements((jintArray)arr2,coldata,0);
//    }
//
//    jint U[n][n];
//    jint L[n][n];
////    jint r[n][n],u[n][n];
////    jint res;
//
////    omp_set_num_threads(n_threads);//设置线程数-----------------------------
//    //计算L、U矩阵
//    for (i = 0; i < n; i++)
//    {
//        U[i][i]=gf_sub(a[i][i],sumLU(i, i, i,(jint **)L,(jint **)U,n));
//        L[i][i] = 1;
////#pragma omp parallel for    //并行化-----------------------------
//        for (j = i+1; j < n; j++)
//        {
//            //按照递推公式进行计算
//            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
//            L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//        }
//
//    }
//    // Ax=b  =>	LUx=b    令Ux=y 有 Ly=b
//    // 1.求出y
//    // 2.已知y  Ux=y  求x
//
////pDatMat 即编码矩阵，每一行进行遍历求解
//    int nDatMatCol;
//    nDatMatCol=col2;
//    jint pIntMat[nDatMatCol][nDatMatCol];
//    jint pResMat[n][n];
////    jint nSum[n];
//
//    for (int nColCnt = 0; nColCnt < nDatMatCol; nColCnt++)
//    {
//        pIntMat[0][nColCnt] = pDatMat[0][nColCnt];
////        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 1 %d %d",n, nDatMatCol);
////        #pragma omp parallel for  //并行化-----------------------------
//        for (i = 1; i<n; i++)
//        {
//
//            for (j = 0; j < n; j++)//
//            {
//                GFType nSum = 0;
////                #pragma omp parallel for  //并行化-----------------------------
//
//                for (k = 0; k < i; k++)
//                {
////                    #pragma omp critical
//                    nSum= gf_add(nSum, gf_mul(L[i][k] , pIntMat[k][nColCnt]));
//                    pIntMat[i][nColCnt] = gf_sub(pDatMat[i][nColCnt] , nSum);
//                }
//            }
//        }
//    }
//    for (int nColCnt = 0; nColCnt < nDatMatCol; nColCnt++)
//    {
//        pResMat[n - 1][nColCnt] = gf_div(pIntMat[n- 1][nColCnt] , U[n - 1][n - 1]);
////        #pragma omp parallel for  //并行化-----------------------------
//        for (i = n - 2; i >= 0; i--)
//        {
//            for (j = 0; j < n; j++)//
//            {
//                GFType nSum = 0;
//                for (k = i + 1; k < n; k++)
//                {
//                    nSum = gf_add(nSum,gf_mul(U[i][k] , pResMat[k][nColCnt]));
//                    pResMat[i][nColCnt] = gf_div(gf_sub(pIntMat[i][nColCnt] , nSum)  , U[i][i]);
//                }
//            }
//        }
//    }
////
//    //返回结果                  //为什么还要另外构造一个新的矩阵？ JAVA用法决定？？？
//    jclass intClass =env->FindClass("[I");
//    jobjectArray mat_inv = env->NewObjectArray(row1,intClass,NULL);
//    for(i=0;i<row1;i++){
//        jintArray intdata = env->NewIntArray(row1);
//        env->SetIntArrayRegion(intdata,0,row1,*&pResMat[i]);
//        env->SetObjectArrayElement(mat_inv,i,intdata);
//        env->DeleteLocalRef(intdata);
//    }
//    return mat_inv;
//}
//////LU串行分解函数！-------------------------------------------------------------


//LU并行分解函数！-------------------------------------------------------------
//矩阵求逆

JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_InverseLU(JNIEnv *env,jclass jclass1,jobjectArray matrix1,jobjectArray matrix2,jint threads){
    int n;
    int k,i,j;
    int flag=1;
    jint s,t;
    //GUI界面设置的线程数
    int n_threads = threads;

    //获取矩阵1的行列数
    int row1 = (int)(env->GetArrayLength(matrix1));//调用函数 得到矩阵长度----行
    jarray arr1 = (jarray)(env->GetObjectArrayElement(matrix1,0));//调用函数 得到矩阵第一行元素构成的矩阵？
    int col1 = (int)(env->GetArrayLength(arr1));//调用函数 得到矩阵长度----列

    //获取矩阵2的行列数
    int row2 = (int)(env->GetArrayLength(matrix2));
    jarray arr2 = (jarray)(env->GetObjectArrayElement(matrix2,0));
    int col2 = (int)(env->GetArrayLength(arr2));

    //判断是否满足可乘条件
    if(col1 != row2){
        return NULL;
    }
    //不是方阵返回null
    if(row1 != col1 || row1 != (int)Java_nc_NcUtils_Rank(env,jclass1,matrix1)){
        return NULL;
    }
    //将输入矩阵转化为函数中矩阵
    n=row1;

    //获取矩阵1
    jint a[row1][col1];
    for(i=0;i<row1;i++){
        arr1 = (jarray)(env->GetObjectArrayElement(matrix1,i));
        jint *coldata = (jint*)(env->GetIntArrayElements((jintArray)arr1,0));
        for(j=0;j<col1;j++){
            a[i][j] = abs(coldata[j]);
        }
        env->ReleaseIntArrayElements((jintArray)arr1,coldata,0);
    }

    //获取矩阵2
    jint pDatMat[row2][col2];
    for(i=0;i<row2;i++){
        arr2 = (jarray)(env->GetObjectArrayElement(matrix2,i));
        jint *coldata = (jint*)(env->GetIntArrayElements((jintArray)arr2,0));
        for(j=0;j<col2;j++){
            pDatMat[i][j] = abs(coldata[j]);
        }
        env->ReleaseIntArrayElements((jintArray)arr2,coldata,0);
    }

    jint U[n][n];
    jint L[n][n];

//    n_threads = 2;
//    omp_set_num_threads(n_threads);//设置线程数-----------------------------

    //计算L、U矩阵------------- 并行LU分解中的编码系数矩阵分k块
    //    该原矩阵被分解为k*k的分块矩阵 不过规模是原来的1/k*k
    int divide=n_threads;
//    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 1 ");

    for(int divideNum=1;divideNum<threads;divideNum++)
    {
        //step1 左上角小块 LU分解     通过公式求解L11和U11
        for (i = n*(divideNum-1)/divide; i < n*divideNum/divide; i++)
        {
            U[i][i]=gf_sub(a[i][i],sumLU(i, i, i,(jint **)L,(jint **)U,n));
            L[i][i] = 1;
            for (j = i+1; j <n*divideNum/divide; j++)  //0~k列
            {
                U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
                L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
            }
        }
//        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 2 ");
        //step2  横向串行 分解U
        for (i = n*(divideNum-1)/divide; i < n*divideNum/divide; i++)
        {
            //这个位置单独计算矩阵U    不可以并行
            for (j = n*divideNum/divide; j < n; j++)   //k~k*k列
            {
                U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
            }
        }
        //step3  纵向并行 分解L
//        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 3 ");
        omp_set_num_threads(divide-divideNum);//设置线程数==bolck数量
//        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 3.1 %d",divide-divideNum);
        #pragma omp parallel for
        for(int block=0;block<divide-divideNum;block++)   //每一小块并行计算
        {
            for (i = n*(divideNum-1)/divide; i < n*divideNum/divide; i++)   //按小块计算
            {
                for (int j = n*(divideNum+block)/divide; j < n*(divideNum+block+1)/divide; j++)
                {
//                    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 3.1 %d %d",a[j][i],U[i][i]);
                    L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
                }
            }
        }
//        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 3.2  ");
    }
//    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 4 ");
//    step4   右下角小块 LU分解
    for (i = n*(divide-1)/divide; i < n; i++)   //3k~4k行
    {
        U[i][i]=gf_sub(a[i][i],sumLU(i, i, i,(jint **)L,(jint **)U,n));
        L[i][i] = 1;
        for (j = i+1; j < n; j++)  //3k~4k列
        {
            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
            L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
        }
    }
//    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 5 ");
    // Ax=b  =>	LUx=b    令Ux=y 有 Ly=b
    // 1.求出y 2.已知y  Ux=y  求x
    //pDatMat 即编码矩阵，每一行进行遍历求解
    int nDatMatCol;
    nDatMatCol=col2;
    jint pIntMat[nDatMatCol][nDatMatCol];
    jint pResMat[n][n];
    for (int nColCnt = 0; nColCnt < nDatMatCol; nColCnt++)
    {
        pIntMat[0][nColCnt] = pDatMat[0][nColCnt];
        // 计算Yi 1-2-3-4-5'
        for(int block=0;block<divide;block++)   //每一小块计算
        {
            for (i = n*block/divide; i<n*(block+1)/divide; i++)
            {
                if(i==0)
                {
                    continue;
                }
                GFType nSum = 0;
                for (k = 0; k < i; k++)
                {
                    nSum= gf_add(nSum, gf_mul(L[i][k] , pIntMat[k][nColCnt]));
                }
                pIntMat[i][nColCnt] = gf_sub(abs(pDatMat[i][nColCnt]) , nSum);
//                __android_log_print(ANDROID_LOG_INFO, "app", "---- in 5.2 i= %d  nSum= %d %d %d",i,nSum,pDatMat[i][nColCnt],pIntMat[i][nColCnt]);
            }
        }
//        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 6 ");
        pResMat[n - 1][nColCnt] = gf_div(pIntMat[n- 1][nColCnt] , U[n - 1][n - 1]);
        // 计算Xi 5-4-3-2-1'
        for(int block=divide-1;block>=0;block--)   //每一小块计算
        {
//            __android_log_print(ANDROID_LOG_INFO, "app", "---- in 6.1 ");
            for (i = n*(block+1)/divide-1; i>=n*block/divide; i--)
            {
//                __android_log_print(ANDROID_LOG_INFO, "app", "---- in 6.15 ");
                if(i==n-1)
                {
                    continue;
                }
//                __android_log_print(ANDROID_LOG_INFO, "app", "---- in 6.2 i= %d ",i);
                GFType nSum = 0;
                for (k = i + 1; k < n; k++)
                {
                    nSum = gf_add(nSum,gf_mul(U[i][k] , pResMat[k][nColCnt]));
                }
//                __android_log_print(ANDROID_LOG_INFO, "app", "---- in 6.3 %d ",U[i][i]);
//                __android_log_print(ANDROID_LOG_INFO, "app", "---- in 6.3 %d ",pIntMat[i][nColCnt]);
                pResMat[i][nColCnt] = gf_div(gf_sub(pIntMat[i][nColCnt] , nSum)  , U[i][i]);
//                __android_log_print(ANDROID_LOG_INFO, "app", "---- in 6.4 ");
            }
        }
    }
//    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 7 ");------------------------------------------
//	// step1: 通过公式求解L11和U11
//    for (i = 0; i < n/divide; i++)   //0~k行
//    {
//        U[i][i]=gf_sub(a[i][i],sumLU(i, i, i,(jint **)L,(jint **)U,n));
//        L[i][i] = 1;
////        n_threads = 8;
//        omp_set_num_threads(n_threads);//设置线程数-----------------------------
//        #pragma omp parallel for//------------------------------------------------------------------
//        for (j = i+1; j <n/divide; j++)  //0~k列
//        {
//            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
//            L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//        }
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 2 ");
//    // step2: 接着计算U12 13 14   这些运算全部在Server#1节点上运行
//    for (i = 0; i < n/divide; i++)   //0~k行
//    {
//        //这个位置单独计算矩阵U   不可以并行
//        for (j = n/divide; j < n; j++)   //k~k*k列
//        {
//            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
//        }
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 3 ");
//    // step3: 接着计算L21 31 41   这些运算分离在不同Server节点上运行
//
//    omp_set_num_threads(3);//设置线程数-----------------------------
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 3.1 n%d divide %d ",n,divide);
//    #pragma omp parallel sections no wait // starts a new team
//    {
//        #pragma omp section    // 计算L21  Server1
//        {
//            for (int i = 0; i < n/divide; i++)
//            {
////                __android_log_print(ANDROID_LOG_INFO, "app", "---- in 3.12n%d divide %d ",n,divide);
//                for (int j = n/divide; j < 2*(n/divide); j++)
//                {
//                    L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
////                    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 3.2 j%d",j);
//                }
//            }
//        }
//        #pragma omp section    // 计算L31  Server2
//        {
//            for (int i = 0; i < n/divide; i++)
//            {
//                for (int j = 2*(n/divide); j < 3*(n/divide); j++)
//                {
//                    L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//                }
//            }
//        }
//        #pragma omp section    // 计算L41  Server3
//        {
//            for (int i = 0; i < n/divide; i++)
//            {
//                for (int j = 3*(n/divide); j < n; j++)
//                {
//                    L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//                }
//            }
//        }
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 4 ");
//    // step4: 通过公式求解L22和U22
////    n_threads = 8;
//    omp_set_num_threads(n_threads);//设置线程数-----------------------------
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 4.1 n%d divide %d ",n,divide);
//    for (i = n/divide; i < 2*(n/divide); i++)   //k~2k行
//    {
////        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 4.10 ");
//        U[i][i]=gf_sub(a[i][i],sumLU(i, i, i,(jint **)L,(jint **)U,n));
////        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 4.11 ");
//        L[i][i] = 1;
//
////        n_threads = 4;    这里不需要设置，和702行左右 重复
////        omp_set_num_threads(n_threads);//设置线程数-----------------------------
//        #pragma omp parallel for
//        for (j = n/divide; j < 2*(n/divide); j++)  //k~2k列
//        {
////            __android_log_print(ANDROID_LOG_INFO, "app", "---- in 4.12 ");
//            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
//            L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//        }
////        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 4.2 ");
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 5 ");
//    // step5: 接着计算U23 24   这些运算全部在Server#1节点上运行
//    for (i = n/divide; i < 2*(n/divide); i++)   //k~2k行
//    {
////        #pragma omp parallel for
//        for (j = 2*(n/divide); j < n; j++)   //2k~k*k列
//        {
//            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
//        }
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 6 ");
//    // step6: 接着计算L32 42   这些运算分离在不同Server节点上运行
//    omp_set_num_threads(2);//设置线程数-----------------------------
//    #pragma omp parallel sections // starts a new team
//    {
//
//        #pragma omp section    // 计算L32  Server1
//        {
//            //__android_log_print(ANDROID_LOG_INFO, "app", "1111111");
//            for (int i = n/divide; i < 2*(n/divide); i++)
//            {
//                for (int j = 2*(n/divide); j < 3*(n/divide); j++)
//                {
//                    L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//                }
//            }
//        }
//        #pragma omp section    // 计算L42  Server2
//        {
//            //__android_log_print(ANDROID_LOG_INFO, "app", "222222222");
//            for (int i = n/divide; i < 2*(n/divide); i++)
//            {
//                for (int j = 3*(n/divide); j < n; j++)
//                {
//                    L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//                }
//            }
//        }
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 7 ");
//    // step7: 通过公式求解L33和U33
////    n_threads = 8;//语句尽量放在外层，执行次数少，运行时间少
//    omp_set_num_threads(n_threads);//设置线程数-----------------------------
//    for (i = 2*(n/divide); i < 3*(n/divide); i++)   //2k~3k行
//    {
//        U[i][i]=gf_sub(a[i][i],sumLU(i, i, i,(jint **)L,(jint **)U,n));
//        L[i][i] = 1;
//        #pragma omp parallel for
//        for (j = 2*(n/divide); j < 3*(n/divide); j++)  //2k~3k列
//        {
//            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
//            L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//        }
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 8 ");
//    // step8: 接着计算U34   这些运算全部在Server#1节点上运行
//    for (i = 2*(n/divide); i < 3*(n/divide); i++)   //2k~3k行
//    {
////        #pragma omp parallel for
//        for (j = 3*(n/divide); j < n; j++)   //3k~k*k列
//        {
//            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
//        }
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 9 ");
//    // step9: 接着计算L43   这些运算分离在不同Server节点上运行
//    for (i = 2*(n/divide); i < 3*(n/divide); i++)  //2k~3k列
//    {
////        #pragma omp parallel for
//        for (j = 3*(n/divide); j < n; j++)//3k~4k行
//        {
//            L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//        }
//    }
////    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 10 ");
//    // step10: 通过公式求解L44和U44
//    for (i = 3*(n/divide); i < n; i++)   //3k~4k行
//    {
//        U[i][i]=gf_sub(a[i][i],sumLU(i, i, i,(jint **)L,(jint **)U,n));
//        L[i][i] = 1;
////        n_threads = 4;  //764行左右 已经设置过了，不重复设置了就
////        omp_set_num_threads(n_threads);//设置线程数-----------------------------
//        #pragma omp parallel for
//        for (j = 3*(n/divide); j < n; j++)  //3k~4k列
//        {
//            U[i][j]=gf_sub(a[i][j],sumLU(i, j, i,(jint **)L,(jint **)U,n));
//            L[j][i]=gf_div(gf_sub(a[j][i],sumLU(j, i, i,(jint **)L,(jint **)U,n)),U[i][i]);
//        }
//    }

//    __android_log_print(ANDROID_LOG_INFO, "app", "---- in 11 ");
    //        __android_log_print(ANDROID_LOG_INFO, "app", "---- in 1 %d %d",n, nDatMatCol);


//    // Ax=b  =>	LUx=b    令Ux=y 有 Ly=b
//    // 1.求出y 2.已知y  Ux=y  求x
//    //pDatMat 即编码矩阵，每一行进行遍历求解
//    int nDatMatCol;
//    nDatMatCol=col2;
//    jint pIntMat[nDatMatCol][nDatMatCol];
//    jint pResMat[n][n];
//
//    omp_set_num_threads(1);//设置线程数-----------------------------
//
//    for (int nColCnt = 0; nColCnt < nDatMatCol; nColCnt++)
//    {
//        pIntMat[0][nColCnt] = pDatMat[0][nColCnt];
//        // 计算Y1'
//        for (i = 1; i<n/divide; i++)  //0~k行
//        {
//            GFType nSum = 0;
//            for (k = 0; k < i; k++)
//            {
//                nSum= gf_add(nSum, gf_mul(L[i][k] , pIntMat[k][nColCnt]));
//            }
//            pIntMat[i][nColCnt] = gf_sub(pDatMat[i][nColCnt] , nSum);
//        }
//        // 计算Y2'
//        for (i = n/divide; i<2*(n/divide); i++)  //k~2k行
//        {
//            GFType nSum = 0;
//            for (k = 0; k < i; k++)
//            {
//                nSum= gf_add(nSum, gf_mul(L[i][k] , pIntMat[k][nColCnt]));
//            }
//            pIntMat[i][nColCnt] = gf_sub(pDatMat[i][nColCnt] , nSum);
//        }
//        // 计算Y3'
//        for (i = 2*(n/divide); i<3*(n/divide); i++)  //2k~3k行
//        {
//            GFType nSum = 0;
//            for (k = 0; k < i; k++)
//            {
//                nSum= gf_add(nSum, gf_mul(L[i][k] , pIntMat[k][nColCnt]));
//            }
//            pIntMat[i][nColCnt] = gf_sub(pDatMat[i][nColCnt] , nSum);
//        }
//        // 计算Y4'
//        for (i = 3*(n/divide); i<n; i++)  //3k~4k行
//        {
//            GFType nSum = 0;
//            for (k = 0; k < i; k++)
//            {
//                nSum= gf_add(nSum, gf_mul(L[i][k] , pIntMat[k][nColCnt]));
//            }
//            pIntMat[i][nColCnt] = gf_sub(pDatMat[i][nColCnt] , nSum);
//        }
//
//        pResMat[n - 1][nColCnt] = gf_div(pIntMat[n- 1][nColCnt] , U[n - 1][n - 1]);
//        // 计算X4
//        for (i = n - 2; i >= 3*(n/divide); i--)  //3k~4k行
//        {
//            GFType nSum = 0;
//            for (k = i + 1; k < n; k++)
//            {
//                nSum = gf_add(nSum,gf_mul(U[i][k] , pResMat[k][nColCnt]));
//            }
//            pResMat[i][nColCnt] = gf_div(gf_sub(pIntMat[i][nColCnt] , nSum)  , U[i][i]);
//        }
//        // 计算X3
//        for (i = 3*(n/divide) - 1; i >= 2*(n/divide); i--)  //2k~3k行
//        {
//            GFType nSum = 0;
//            for (k = i + 1; k < n; k++)
//            {
//                nSum = gf_add(nSum,gf_mul(U[i][k] , pResMat[k][nColCnt]));
//            }
//            pResMat[i][nColCnt] = gf_div(gf_sub(pIntMat[i][nColCnt] , nSum)  , U[i][i]);
//        }
//        // 计算X2
//        for (i = 2*(n/divide) - 1; i >= (n/divide); i--)  //k~2k行
//        {
//            GFType nSum = 0;
//            for (k = i + 1; k < n; k++)
//            {
//                nSum = gf_add(nSum,gf_mul(U[i][k] , pResMat[k][nColCnt]));
//            }
//            pResMat[i][nColCnt] = gf_div(gf_sub(pIntMat[i][nColCnt] , nSum)  , U[i][i]);
//        }
//        // 计算X1
//        for (i = (n/divide) - 1; i >= 0; i--)  //0~k行
//        {
//            GFType nSum = 0;
//            for (k = i + 1; k < n; k++)
//            {
//                nSum = gf_add(nSum,gf_mul(U[i][k] , pResMat[k][nColCnt]));
//            }
//            pResMat[i][nColCnt] = gf_div(gf_sub(pIntMat[i][nColCnt] , nSum)  , U[i][i]);
//        }
//    }
    //返回结果                  //为什么还要另外构造一个新的矩阵？ JAVA用法决定？？？
    jclass intClass =env->FindClass("[I");
    jobjectArray mat_inv = env->NewObjectArray(row1,intClass,NULL);
    for(i=0;i<row1;i++){
        jintArray intdata = env->NewIntArray(row1);
        env->SetIntArrayRegion(intdata,0,row1,*&pResMat[i]);
        env->SetObjectArrayElement(mat_inv,i,intdata);
        env->DeleteLocalRef(intdata);
    }
    return mat_inv;
}
////LU并行分解函数！-------------------------------------------------------------
//



//有限域求方阵的行列式的值
JNIEXPORT jint JNICALL
Java_nc_NcUtils_Det(JNIEnv *env,jclass ,jobjectArray matrix){
    int i,j,z,x;

    //获取矩阵1的行列数
    int row = (int)(env->GetArrayLength(matrix));
    jarray arr = (jarray)(env->GetObjectArrayElement(matrix,0));
    int col = (int)(env->GetArrayLength(arr));
    env->DeleteLocalRef(arr);

    //判断是否是方阵
    if(col != row){
        return NULL;
    }

    //获取矩阵
    jint mat[row][col];
    for(i=0;i<row;i++){
        arr = (jarray)(env->GetObjectArrayElement(matrix,i));
        jint *data = (jint*)(env->GetIntArrayElements((jintArray)arr,0));
        for(j=0;j<col;j++){
            mat[i][j] = data[j];
        }
        env->ReleaseIntArrayElements((jintArray)arr,data,0);
    }



    //将矩阵转换成下三角矩阵
    for (i = 0; i < row; i++) {
        if (mat[i][i] == 0) {
            for (z = i + 1; z < row; z++) {
                if (mat[z][i] != 0) {
                    for (x = 0; x < row; x++) {
                        int temp = mat[i][x];
                        mat[i][x] = mat[z][x];
                        mat[z][x] = temp;
                    }
                    break;
                }
            }
        }
        for (j = i + 1; j < row; j++) {
            GFType temp = gf_div(mat[j][i], mat[i][i]);
            for (z = 0; z < row; z++) {
                mat[j][z] = gf_add(mat[j][z], gf_mul(temp, mat[i][z]));
            }
        }
    }

    //计算行列式的值
    jint determinant = mat[0][0];
    for(i=1;i<row;i++){
        determinant = gf_mul(determinant,mat[i][i]);
    }
    return determinant;
}

JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_SparseMatrix(JNIEnv *env,jclass ,jint n,jint row,jint col){
    //创建一个新矩阵 长度为row*col
    jint ** mat;
    int i,j,temp,flag;
    mat = new jint *[row];
    for(i=0;i<row;i++){
        mat[i] = new jint[col];
    }
    //将矩阵赋值为零矩阵
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            mat[i][j] = 0;
        }
    }
    //通过随机函数和稀疏度进行赋值
    for(i=0;i<row;i++){
        for(j=0;j<n;j++){
            flag = 1;
            while(flag){
                temp = rand()%col;
                if(mat[i][temp] == 0){
                    mat[i][temp] = (jint)gf();
                    flag = 0;
                }
            }

        }
    }
    //将生成的矩阵传回java层
    jclass intClass =env->FindClass("[I");
    jobjectArray matrix = env->NewObjectArray(row,intClass,NULL);
    for(i=0;i<row;i++){
        jintArray intdata = env->NewIntArray(col);
        env->SetIntArrayRegion(intdata,0,col,*&mat[i]);
        env->SetObjectArrayElement(matrix,i,intdata);
        env->DeleteLocalRef(intdata);
    }
    return matrix;
}
/*JNIEXPORT jobjectArray JNICALL
Java_nc_NcUtils_SparseMatrix1(JNIEnv *env,jclass ,jint n,jint col){

}*/




