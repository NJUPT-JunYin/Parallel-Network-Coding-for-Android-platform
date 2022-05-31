package com.hl.testnc;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;

import com.example.testnc.R;

import nc.NcUtils;

public class MatrixMultiplyActivity extends AppCompatActivity {

    private Button btn_CreateMatrix;
    private Button btn_Multiply;
    private int[][] matrix1;
    private int[][] matrix2;
    private int[][] matrix3;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_matrix_multiply);
        NcUtils.initGalois(8);          //界面生成的时候进行有限域初始化

        btn_CreateMatrix = findViewById(R.id.btn_generate);
        btn_Multiply = findViewById(R.id.btn_Multiply);

        btn_CreateMatrix.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
//                matrix1 = MatrixMultiplyActivity.generateMatrix(3,4);
                matrix1 = NcUtils.Matrix(5,5);
//                matrix1 = NcUtils.Matrix(3,3);
//                matrix2 = MatrixMultiplyActivity.generateMatrix(4,3);
                matrix2 = NcUtils.Matrix(5,10);

                TextView tv_matrix1 = findViewById(R.id.tv_matrix1);
                MatrixMultiplyActivity.printMatrix(tv_matrix1,matrix1);

                TextView tv_matrix2 = findViewById(R.id.tv_matrix2);
                MatrixMultiplyActivity.printMatrix(tv_matrix2,matrix2);
            }
        });

        btn_Multiply.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                matrix3 = NcUtils.MatrixMultiply(matrix1,matrix2);
//                matrix3 = NcUtils.SparseMatrix(3,6,6);
//                matrix3 = NcUtils.Matrix(3,3);
                TextView tv_matrix3 = findViewById(R.id.tv_matrix3);
                MatrixMultiplyActivity.printMatrix(tv_matrix3,matrix3);
//
//                MatrixMultiplyActivity.printMatrix(tv_matrix3,NcUtils.Inverse(matrix1));
//                tv_matrix3.setText(Integer.toString(NcUtils.gf()));
//                tv_matrix3.setText(Integer.toString(NcUtils.Det(matrix1)));
            }
        });
    }

    @Override
    protected void onDestroy() {
        super.onDestroy();
        NcUtils.uninitGalois();     //界面退出的时候执行反初始化操作避免内存泄漏
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
