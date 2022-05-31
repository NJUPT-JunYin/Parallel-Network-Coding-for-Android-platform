package com.hl.testnc;

import androidx.appcompat.app.AlertDialog;
import androidx.appcompat.app.AppCompatActivity;

import android.content.Context;

import android.os.Bundle;
import android.view.MotionEvent;
import android.view.View;
import android.view.inputmethod.InputMethodManager;
import android.widget.AdapterView;
import android.widget.Button;
import android.widget.EditText;

import android.widget.ListView;
import android.widget.SimpleAdapter;
import android.widget.TextView;
import android.widget.Toast;

import com.example.testnc.R;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.text.SimpleDateFormat;

import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.lang.Math;


import nc.NcUtils;

import static java.lang.System.nanoTime;

public class NetworkCodingActivity extends AppCompatActivity {

    private Context mContext;
    private Button btn_InitGf;
    private Button btn_GenerateFile;
    private Button btn_OpenFile;
    private EditText tv_GfM;
    private EditText et_FileSize;
    private EditText et_sparsity;//输入稀疏度
    private String filename = "";
    private byte[] file;
    private Button btn_coding;
    private Button btn_decoding;
    private Button btn_sparsity;//稀疏度
    private int gfsize;
    private int psize;
    private int block;
    private int countmax;
    //private int[][] matrix;
    //private long time1;
    //private long time2;
    private long sum=0;
    private long EncodingTime;
    //private  int getblocknum=0;
    //private int count=0;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_network_coding);
        mContext = getApplicationContext();

        //初始化有限域的代码
        btn_InitGf = findViewById(R.id.btn_initgf);
        btn_InitGf.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                tv_GfM = findViewById(R.id.et_gfm);
                int m = Integer.parseInt(tv_GfM.getText().toString());
                NcUtils.initGalois(m);
                gfsize = 2<<(m-1);
                Toast toast = Toast.makeText(mContext,"Initialization of finite field GF（2^"+ m +")succeeded",Toast.LENGTH_SHORT);
                toast.show();
            }
        });

        //生成随机大小文件的代码
        btn_GenerateFile = findViewById(R.id.btn_generatefile);
        btn_GenerateFile.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                et_FileSize = findViewById(R.id.et_filesize);
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

                Toast toast = Toast.makeText(mContext,"File generation succeeded，File size is"+file.length+"B",Toast.LENGTH_SHORT);
                toast.show();
            }
        });

        //打开某个文件的代码
        btn_OpenFile = findViewById(R.id.btn_openfile);
        btn_OpenFile.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                TextView tv_fileStatue = findViewById(R.id.tv_fileStatue);
                tv_fileStatue.setText("Currently selected file:"+filename+",\nFile size："+file.length+"B");
            }
        });

        //选择稀疏度的代码
        btn_sparsity = findViewById(R.id.btn_sparsity);
        btn_sparsity.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                et_sparsity = findViewById(R.id.et_sparsity);
                int n=Integer.parseInt((et_sparsity.getText().toString()));
                countmax = n*20/100;
                Toast toast = Toast.makeText(mContext,"Sparsity is selected as "+ countmax*100/20 +"%successful",Toast.LENGTH_SHORT);
                toast.show();
            }
        });
        //网络编码按钮
        btn_coding =  findViewById(R.id.btn_networkcoding);
        btn_coding.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(filename == ""){
                    return;
                }

               //获取设置的文件块大小和分块数
                EditText et_blocklen = findViewById(R.id.et_blocksize);
                psize = Integer.parseInt(et_blocklen.getText().toString());//文件块长度
                block = file.length / psize;//编码块数量=文件长度/文件块长度

                //待编码文件读取至内存
                //long time1 = System.nanoTime();//开始计时
                int[][] codingblock = new int[countmax][psize/4];
                int matrix1[]=new int[block];
                for(int runtime=0;runtime<1000;runtime++){
                for(int r=0;r<block;r++){
                       matrix1[r]=0;
                  }
                 for(int r=0;r<countmax;r++){
                       int temp=(int)(Math.random()*block+0);
                       matrix1[temp]=(int)(Math.random()*255+1);
                 }
                        int getblocknum=0;
                        long time1 = System.nanoTime();//开始计时
                        //Toast toast = Toast.makeText(mContext,"编码时间为"+ time1 ,Toast.LENGTH_SHORT);
                        //toast.show();
                       for(int i=0;i<block;i++){
                            if(matrix1[i]==0){
                                continue ;
                            }
                            for(int j=0;j<psize/4;j++){
                                //getblocknum++;
                                codingblock[getblocknum][j] = ((int)((file[j*4] & 0xFF)//1int=4B;1B=8bit;
                                        | ((file[1+j*4] & 0xFF)<<8)
                                        | ((file[2+j*4] & 0xFF)<<16)
                                        | ((file[3+j*4] & 0xFF)<<24))) % gfsize;
                            }
                            getblocknum++;
                            if(getblocknum==countmax){
                                    break ;
                            }
                        }
                        int count = 0;      //已编码块数
                        int[] rowVector = new int[psize/4];//行向量存int型；1（int）=4B，int占4个字节
                        for(int j = 0; j < block;j++ ){
                            int codingnum = matrix1[j];//改成向量
                            if(codingnum == 0)
                                continue;
                            count++;
                            for(int k = 0; k < psize/4 ; k++){
                                rowVector[k] += codingblock[count-1][k] * codingnum;
                            }
                            if(count == countmax)//稀疏度
                                break ;
                        }
                long time2 = System.nanoTime();
                long EncodingTime = time2 - time1;
                sum+=EncodingTime;
                }
                sum/=1000;
                //TextView tv_result = findViewById(R.id.tv_result);
                //tv_result.setText("编码运算时间："+Long.toString(sum)+"ns");
            }
        });

        //解码测试按钮
        btn_decoding = findViewById(R.id.btn_decoding);
        btn_decoding.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                EditText et_blocklen = findViewById(R.id.et_blocksize);
                psize = Integer.parseInt(et_blocklen.getText().toString());//文件块长度
                block = file.length / psize;//编码块数量=文件长度/文件块长度//block=20;
                int i,j,k,q,m;
                float [][]matrix2=new float[block][psize/4];
                int [][]matrix=new int[block][psize/4];
                do{
                    matrix = NcUtils.SparseMatrix(countmax,block,block);
                }while(NcUtils.Rank(matrix) != block);//生成满秩稀疏矩阵*/
                for(i=0;i<20;i++){
                    for(j=0;j<20;j++){
                        matrix2[i][j]=(float)matrix[i][j];
                    }
                }

                int [][]codingfile=NcUtils.Matrix(20,psize/4);//生成随机已编码文件；

                //float [][]Inversematrix=NcUtils.LUInverse(matrix2,20);//求稀疏矩阵的逆；
                int N=block;
                //int i,j,k;

                float [][]L=new float[N][N];
                float [][]U=new float[N][N];
                float [][]out=new float[N][N];
                float [][]r=new float[N][N];
                float [][]u=new float[N][N];
                //out矩阵保存矩阵a的逆
                //r为L矩阵的逆u为U矩阵的逆
                int n=N;
                int flag=1;
                float s,t;

                for(j=0;j<n;j++)
                    matrix2[0][j]=matrix2[0][j];  //计算U矩阵的第一行

                for(i=1;i<n;i++)
                    matrix2[i][0]=matrix2[i][0]/matrix2[0][0];   //计算L矩阵的第1列
                for(k=1;k<n;k++)
                {
                    for(j=k;j<n;j++)
                    {
                        s=0;
                        for (i=0;i<k;i++)
                            s=s+matrix2[k][i]*matrix2[i][j];   //累加
                        matrix2[k][j]= (matrix2[k][j]-s); //计算U矩阵的其他元素
                    }
                    for(i=k+1;i<n;i++)
                    {
                        t=0;
                        for(j=0;j<k;j++)
                            t=t+matrix2[i][j]*matrix2[j][k];   //累加
                        matrix2[i][k]=(matrix2[i][k]-t)/matrix2[k][k];    //计算L矩阵的其他元素
                    }
                }
                for(i=0;i<n;i++)
                    for(j=0;j<n;j++)
                    {
                        if(i>j)
                        {
                            L[i][j]=matrix2[i][j];
                            U[i][j]=0;
                        }//如果i>j，说明行大于列，计算矩阵的下三角部分，得出L的值，U的//为0
                        else
                        {
                            U[i][j]=matrix2[i][j];
                            if(i==j)
                                L[i][j]=1;  //否则如果i<j，说明行小于列，计算矩阵的上三角部分，得出U的值，L的为0
                            else
                                L[i][j]=0;
                        }
                    }
                if(U[1][1]*U[2][2]*U[3][3]*U[4][4]==0) {
                    flag = 0;
                    //printf("\n逆矩阵不存在");}
                }
                if(flag==1)
                {
                    for (i=0;i<n;i++) /*求矩阵U的逆 */
                    {
                        u[i][i]=1/U[i][i];//对角元素的值，直接取倒数
                        for (k=i-1;k>=0;k--)
                        {
                            s=0;
                            for (j=k+1;j<=i;j++)
                                s=s+U[k][j]*u[j][i];
                            u[k][i]=-s/U[k][k];//迭代计算，按列倒序依次得到每一个值，
                        }
                    }
                    for (i=0;i<n;i++) //求矩阵L的逆
                    {
                        r[i][i]=1; //对角元素的值，直接取倒数，这里为1
                        for (k=i+1;k<n;k++)
                        {
                            for (j=i;j<=k-1;j++)
                                r[k][i]=r[k][i]-L[k][j]*r[j][i];   //迭代计算，按列顺序依次得到每一个值
                        }
                    }

                    //////////将r和u相乘，得到逆矩阵
                    for(i=0;i<n;i++)
                    {
                        for(j=0;j<n;j++)
                        {out[i][j]=0;}
                    }
                    for(i=0;i<n;i++)
                    {
                        for(j=0;j<n;j++)
                        {
                            for(k=0;k<n;k++)
                            {
                                out[i][j]+=u[i][k]*r[k][j];
                            }
                        }
                    }
                }

                float [][]decodingfile = new float[block][psize/4];

                float []sum=new float[psize/4];

                long time1=System.nanoTime();
                for(i=0;i<block;i++){
                    for(m=0;m<psize/4;m++){
                        sum[m]=0;
                    }
                    for(j=0;j<block;j++){
                        if(/*Inversematrix*/out[i][j]==0)
                            continue;
                        else {
                            float codingnum=/*Inversematrix*/out[i][j];
                            for(k=0;k<psize/4;k++){
                                sum[k]+=codingnum*codingfile[j][k];
                            }
                        }
                    }
                    for(q=0;q<psize/4;q++){
                        decodingfile[i][q]=sum[q];
                    }
                }
                long time2=System.nanoTime();
                long time0=time2-time1;
                TextView tv_result = findViewById(R.id.tv_result);
                tv_result.setText("Decoding operation time："+Long.toString(time0)+"ns");
            }
        });
    }

    //弹出对话框
    private void showDialog(){
        AlertDialog.Builder builder = new AlertDialog.Builder(this);
        String[] files = mContext.fileList();
        ListView listView = new ListView(this);
        listView.setFadingEdgeLength(0);
        List<Map<String,String>> list = new ArrayList<>();
        for(int m=0;m<files.length;m++){
            Map<String,String> map = new HashMap<>();
            map.put("filename",files[m]);
            list.add(map);
        }
        SimpleAdapter adapter = new SimpleAdapter(mContext, list, R.layout.openfileitem,
                new String[]{"filename"},new int[]{R.id.tv_item});
        listView.setAdapter(adapter);
        builder.setView(listView);
        builder.setTitle("选择文件");
        final AlertDialog dialog = builder.create();
        dialog.show();
        listView.setOnItemClickListener(new AdapterView.OnItemClickListener() {
            @Override
            public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
                TextView tv = (TextView) view.findViewById(R.id.tv_item);
                filename = tv.getText().toString();
                TextView tv_filestatue = findViewById(R.id.tv_fileStatue);
                tv_filestatue.setText("当前选择文件："+filename);
                dialog.cancel();
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


}
