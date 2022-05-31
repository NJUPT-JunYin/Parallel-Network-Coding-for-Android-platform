package com.hl.testnc;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;

import com.example.testnc.R;


public class MainActivity extends AppCompatActivity {

    /**
     * 1kb - 4kb 文件块
     * @param savedInstanceState
     */

    private Button btn_MatrixMultiply;
    private Button btn_Coding;
    private Button btn_Coding2;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        btn_MatrixMultiply = findViewById(R.id.btn_MartixMutiply);
        btn_MatrixMultiply.setOnClickListener(new View.OnClickListener(){

            @Override
            public void onClick(View view) {
                Intent intent = new Intent(MainActivity.this,MatrixMultiplyActivity.class);
                startActivity(intent);
            }
        });

        btn_Coding = findViewById(R.id.btn_coding);
        btn_Coding.setOnClickListener(new View.OnClickListener(){

            @Override
            public void onClick(View view) {
                Intent intent = new Intent(MainActivity.this,NetworkCodingActivity.class);
                startActivity(intent);
            }
        });

        btn_Coding2 = findViewById(R.id.btn_coding2);
        btn_Coding2.setOnClickListener(new View.OnClickListener(){

            @Override
            public void onClick(View view) {
                Intent intent = new Intent(MainActivity.this,NetworkCoding2Activity.class);
                startActivity(intent);
            }
        });
    }

}
