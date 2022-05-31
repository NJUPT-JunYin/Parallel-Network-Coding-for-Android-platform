package com.hl.testnc;

import androidx.appcompat.app.AppCompatActivity;
import android.content.Context;
import android.os.Bundle;
import android.view.MotionEvent;
import android.view.View;
import android.view.inputmethod.InputMethodManager;

import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;
import android.widget.Toast;

import com.example.testnc.R;

import java.text.SimpleDateFormat;

import java.util.Date;

import nc.NcUtils;

public class NetworkCoding2Activity extends AppCompatActivity {

    private Context mContext;
    private Button btn_InitGf;
    private Button btn_GenerateFile;
    private Button btn_OpenFile;
    private EditText tv_GfM;
    private EditText et_FileSize;
    private String filename = "";
    private Button btn_coding;
    private Button btn_decoding;
    private byte[] file;
    private int gfsize;
    private int[][] matrix;
    private int[][] result;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_network_coding2);
        mContext = getApplicationContext();

        //初始化有限域的代码
        btn_InitGf = findViewById(R.id.btn_initgf2);
        btn_InitGf.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                tv_GfM = findViewById(R.id.et_gfm2);
                int m = Integer.parseInt(tv_GfM.getText().toString());
                NcUtils.initGalois(8);
                gfsize = 2<<(m-1);
                    Toast toast = Toast.makeText(mContext,"Initial GF（"+gfsize+")",Toast.LENGTH_SHORT);
                toast.show();
            }
        });

        //生成随机大小文件的代码
        btn_GenerateFile = findViewById(R.id.btn_generatefile2);
        btn_GenerateFile.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                et_FileSize = findViewById(R.id.et_filesize2);
                int size = Integer.parseInt((et_FileSize.getText().toString()));

                //根据系统时间生成文件名
                Date date = new Date();
                SimpleDateFormat sdf = new SimpleDateFormat("HHmmss");
                filename = "F" + sdf.format(date) + ".txt";

                int filelen =0;
                file = new byte[size];
                //新建一个文件并随机生成ascii写入文件
                while(filelen<size){
                    file[filelen++] = (byte)(char)(int)(Math.random() * 10000 % 256);
                }

                Toast toast = Toast.makeText(mContext,"Sucess，size"+file.length+"B",Toast.LENGTH_SHORT);
                toast.show();
            }
        });

//        //打开某个文件的代码
//        btn_OpenFile = findViewById(R.id.btn_openfile2);
//        btn_OpenFile.setOnClickListener(new View.OnClickListener() {
//            @Override
//            public void onClick(View view) {
//                TextView tv_fileStatue2 = findViewById(R.id.tv_fileStatue2);
//                tv_fileStatue2.setText("当前选择的文件:"+filename+",\n文件大小："+file.length+"B");
//            }
//        });


        //网络编码按钮
        btn_coding =  findViewById(R.id.btn_networkcoding2);
        btn_coding.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(filename.equals("")){
                    return;
                }

                //获取设置的文件块大小和分块数
                EditText et_blocklen = findViewById(R.id.et_blocksize2);
                int psize = Integer.parseInt(et_blocklen.getText().toString());
                int block = file.length / psize;

                //获取结果输出框
                TextView tv_result = findViewById(R.id.tv_result2);

                long startTime = System.currentTimeMillis();
                int[][] codingblock = new int[block][psize/4];


                for(int i=0;i<block;i++){
                    for(int j=0;j<psize/4;j++){
                        codingblock[i][j] = ((int)((file[j*4] & 0xFF)
                                | ((file[1+j*4] & 0xFF)<<8)
                                | ((file[2+j*4] & 0xFF)<<16)
                                | ((file[3+j*4] & 0xFF)<<24))) % gfsize;
                    }
                }
//                int[][] codingblock = new int[block][psize];
//                for(int i=0;i<block;i++){
//                    for(int j=0;j<psize;j++){
//                        codingblock[i][j] = file[i*psize +j];
//                    }
//                }
                //生成有限域上的随机矩阵
                long time1 = System.currentTimeMillis();
                matrix = NcUtils.Matrix(block,block);
                long time2 = System.currentTimeMillis();
                double time0 = time2 - time1;

                //获取编码块
                result = NcUtils.MatrixMultiply(matrix,codingblock);

                long endTime = System.currentTimeMillis();

                long runTime = endTime - startTime;
                double avgTime = ((double)runTime / block) * 1000;


                tv_result.setText("Total time of coding："+Long.toString(runTime)+"ms\n" +
                        "Coding time of generating a coded bolck："+Double.toString(avgTime)+"us\n" +
                        "Time of generating the coding martix："+Double.toString(time0)+"ms");
            }
        });

        //解码测试按钮
        btn_decoding = findViewById(R.id.btn_decoding2);
        btn_decoding.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                //求逆矩阵 然后左乘解码

                int psize = result[0].length*4;
                int block = result.length;

//                int[][] decode = new int[block][psize];
//                for(int i=0;i<block;i++){
//                    for(int j=0;j<psize/4;j++){
//                        int temp = result[i][j];
//                        decode[i][j*4] = (int)(byte)(temp & 0xFF);
//                        decode[i][j*4+1] = (int)(byte)((temp>>8) & 0xFF);
//                        decode[i][j*4+2] = (int)(byte)((temp>>16) & 0xFF);
//                        decode[i][j*4+3] = (int)(byte)((temp>>24) & 0xFF);
//                    }
//                }
                int[][] decode = new int[block][psize/4];
                for(int i=0;i<block;i++){
                    for(int j=0;j<psize/4;j++) {
                        decode[i][j] = result[i][j] % gfsize;
                    }
                }

                TextView tv_result = findViewById(R.id.tv_result4);


                int[][] inverse=null;
                long startTime = System.currentTimeMillis();

                et_FileSize = findViewById(R.id.et_filesize5);
                int threads = Integer.parseInt((et_FileSize.getText().toString()));

//                int[][] decoding1 = NcUtils.InverseLU(matrix, decode);  //LU矩阵解码-----------------------
//                for(int i=0;i<1000;i++){
//                    int[][] decoding1 = NcUtils.InverseLU(matrix, decode);
//
                for(int i=0;i<1;i++){
                    int[][] decoding1 = NcUtils.InverseLU(matrix, decode,threads);
                }
                long invTimeLU = System.currentTimeMillis();
                double runTimeLU = (double)(invTimeLU - startTime)/1;

//                tv_result.setText("解码中请等待--- \n");
//                inverse = NcUtils.Inverse(matrix);  //高斯矩阵求逆①-----------------------
//                for(int i=0;i<1000;i++){
//                    inverse= NcUtils.Inverse(matrix);
//                }
//高斯need------------------------------------------
                for(int i=0;i<1;i++){
                    inverse= NcUtils.Inverse(matrix);
                }
                long invTime = System.currentTimeMillis();
                double runTime1 = (double)(invTime - invTimeLU)/1;
//高斯need------------------------------------------

//                for(int i=0;i<1000;i++){
//                    int[][] decoding2 = NcUtils.MatrixMultiply(inverse , decode);//解码 ②-----------------------
//                }

//高斯need------------------------------------------
                for(int i=0;i<1;i++){
                    int[][] decoding2 = NcUtils.MatrixMultiply(inverse , decode);//解码 ②-----------------------
                }
//高斯need------------------------------------------
                long endTime = System.currentTimeMillis();
                double runTime2 = (double)(endTime - invTime)/1;


//
//
//                double avgTime = ((double)(runTime1 /1000 + runTime2) / block) * 1000;
//                double avgTime = ((double)(runTime1  + runTime2) / block) * 1000;
                TextView tv_result2 = findViewById(R.id.tv_result3);
                tv_result2.setText("Decode ending!    Threads: "+(threads)+"\n"+
//                                (threads)+"\n" +
//                        "平均每个文件块解码时间："+Double.toString(avgTime)+"us\n"+
//                        "1000次平均 高斯求逆时间:"+Double.toString(runTime1)+"ms\n"+
//                        "1000次平均 高斯解码总时间："+Double.toString(runTime1+ runTime2)+"ms\n" +
//                        "1000次平均 多线程并行LU解码时间:"+Double.toString(runTimeLU)+"ms\n");
//                        "高斯求逆时间:"+Double.toString(runTime1)+"ms\n"+
                        "Time of gauss decoding："+Double.toString(runTime1+ runTime2)+"ms\n" +
                        "Time of parallel LU decoding:"+Double.toString(runTimeLU)+"ms\n");

            }
        });

    }

    //用于关闭输入法
    @Override
    public boolean dispatchTouchEvent(MotionEvent ev) {
        if (ev.getAction() == MotionEvent.ACTION_DOWN) {
            View v = getCurrentFocus();
            if (isShouldHideInput(v, ev)) {

                InputMethodManager imm = (InputMethodManager) getSystemService(Context.INPUT_METHOD_SERVICE);
                if (imm != null) {
                    imm.hideSoftInputFromWindow(v.getWindowToken(), 0);
                }
            }
            return super.dispatchTouchEvent(ev);
        }
        if (getWindow().superDispatchTouchEvent(ev)) {
            return true;
        }
        return onTouchEvent(ev);
    }

    private boolean isShouldHideInput(View v, MotionEvent event) {
        if (v != null && (v instanceof EditText)) {
            int[] location = {0, 0};
            v.getLocationOnScreen(location);
            int left = location[0];
            int top = location[1];

            if (event.getX() < left || (event.getX() > left + v.getWidth())
                    || event.getY() < top || (event.getY() > top + v.getHeight())) {
                return true;
            } else {
                return false;
            }
        }
        return false;
    }

    //关闭界面时反初始化有限域
    @Override
    protected void onDestroy() {
        super.onDestroy();
        NcUtils.uninitGalois();
    }


    public static void printMatrix(TextView tv, int[][] matrix){
        String str = "";
        for(int i=0;i<matrix.length;i++){
            for(int j=0;j<matrix[i].length;j++){
                str += matrix[i][j] + " ";
            }
            str += "\n";
        }
        tv.setText(str);
    }
}

