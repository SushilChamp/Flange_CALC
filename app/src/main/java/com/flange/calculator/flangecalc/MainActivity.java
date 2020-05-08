package com.flange.calculator.flangecalc;

import android.content.pm.ResolveInfo;
import android.graphics.Paint;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.CompoundButton;
import android.widget.EditText;
import android.widget.TextView;
import android.widget.Toast;
import android.widget.ToggleButton;
import android.Manifest;
import android.content.DialogInterface;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.net.Uri;
import android.os.Build;
import android.os.Environment;
import android.support.v4.app.ActivityCompat;
import android.support.v7.app.AlertDialog;
import android.util.Log;

import com.flange.calculator.flangecalc.R;
import com.itextpdf.text.BaseColor;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Element;
import com.itextpdf.text.Paragraph;
import com.itextpdf.text.pdf.PdfPCell;
import com.itextpdf.text.pdf.PdfPTable;
import com.itextpdf.text.pdf.PdfWriter;
import com.itextpdf.text.Font;
import com.itextpdf.text.PageSize;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.Image;
import com.itextpdf.text.Phrase;
import com.itextpdf.text.html.WebColors;
import com.itextpdf.text.ListItem;

import java.io.DataInput;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.sql.Timestamp;
import java.util.Date;
import java.text.DateFormat;
import java.time.format.TextStyle;
import java.lang.Math;
import java.text.SimpleDateFormat;
import java.util.List;
import java.util.ArrayList;
import java.math.RoundingMode;
import java.text.DecimalFormat;

//this is starting of main function
public class MainActivity extends AppCompatActivity {

    //We declare global variables in this section
    EditText pressure;
    EditText temp;
    EditText force;
    EditText moment;
    EditText gasketod;
    EditText gasketid;
    EditText gasketcontactarea;
    EditText flangeod;
    EditText flangeid;
    EditText flangebcd;
    EditText flangethk;
    EditText hubsmall;
    EditText hublarge;
    EditText hublength;
    EditText rigidity;
    EditText boltnum;
    EditText boltdia;
    Button submit;
    ToggleButton ptype;
    ToggleButton gtype;
    ToggleButton ftype;
    //TextView resultN;
    private File pdfFile;
    final private int REQUEST_CODE_ASK_PERMISSIONS = 111;


    //Function #1----assigns values to the global variables, takes toggle buttons, check for empty and calls function for calculation
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);


        //Underlining texts in front end
        TextView designdata = findViewById(R.id.designdata);
        designdata.setPaintFlags(designdata.getPaintFlags()| Paint.UNDERLINE_TEXT_FLAG);
        TextView gasketdata = findViewById(R.id.gasketdata);
        gasketdata.setPaintFlags(gasketdata.getPaintFlags()| Paint.UNDERLINE_TEXT_FLAG);
        TextView flangedata = findViewById(R.id.flangedata);
        flangedata.setPaintFlags(flangedata.getPaintFlags()| Paint.UNDERLINE_TEXT_FLAG);


        //we assign values from user input to the global variables here
        pressure = (EditText) findViewById(R.id.pressure);
        temp = (EditText) findViewById(R.id.temp);
        force = (EditText) findViewById(R.id.force);
        moment = (EditText) findViewById(R.id.moment);
        gasketod = (EditText) findViewById(R.id.gasketod);
        gasketid = (EditText) findViewById(R.id.gasketid);
        gasketcontactarea = (EditText) findViewById(R.id.gasketcontactarea);
        flangeod = (EditText) findViewById(R.id.flangeod);
        flangeid = (EditText) findViewById(R.id.flangeid);
        flangebcd = (EditText) findViewById(R.id.flangebcd);
        flangethk = (EditText) findViewById(R.id.flangethk);
        hubsmall = (EditText) findViewById(R.id.hubsmall);
        hublarge = (EditText) findViewById(R.id.hublarge);
        hublength = (EditText) findViewById(R.id.hublength);
        rigidity = (EditText) findViewById(R.id.rigidity);
        boltnum = (EditText) findViewById(R.id.boltnum);
        boltdia = (EditText) findViewById(R.id.boltdia);
        submit = (Button) findViewById(R.id.butsubmit);
        ptype = (ToggleButton) findViewById(R.id.butpressure);
        gtype = (ToggleButton) findViewById(R.id.butgasket);
        ftype = (ToggleButton) findViewById(R.id.butflange);
        //resultN = (TextView) findViewById(R.id.resultN);




        ptype.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                if(isChecked) {
                    Toast.makeText(getApplicationContext(), "External is selected!", Toast.LENGTH_LONG).show();
                }else{
                    Toast.makeText(getApplicationContext(),"Internal is selected!",Toast.LENGTH_LONG).show(); } }});




        gtype.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                if(isChecked) {
                    Toast.makeText(getApplicationContext(), "SEG is selected!", Toast.LENGTH_LONG).show();
                }else{
                    Toast.makeText(getApplicationContext(),"NSEG is selected!",Toast.LENGTH_LONG).show(); } }});




        ftype.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                if(isChecked) {
                    Toast.makeText(getApplicationContext(), "Loose is selected!", Toast.LENGTH_LONG).show();
                }else{
                    Toast.makeText(getApplicationContext(),"Integral is selected!",Toast.LENGTH_LONG).show(); } }});



        //This is where we check for impty inputs and call other function "calculate flange" for further calculation
        submit.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {

                if (boltdia.getText().toString().isEmpty()) {
                    boltdia.setError("Fill this");
                    boltdia.requestFocus();}
                if (boltnum.getText().toString().isEmpty()) {
                    boltnum.setError("Fill this");
                    boltnum.requestFocus();}
                if (hublength.getText().toString().isEmpty()) {
                    hublength.setError("Fill this");
                    hublength.requestFocus();}
                if (hublarge.getText().toString().isEmpty()) {
                    hublarge.setError("Fill this");
                    hublarge.requestFocus();}
                if (hubsmall.getText().toString().isEmpty()) {
                    hubsmall.setError("Fill this");
                    hubsmall.requestFocus();}
                if (flangethk.getText().toString().isEmpty()) {
                    flangethk.setError("Fill this");
                    flangethk.requestFocus();}
                if (flangebcd.getText().toString().isEmpty()) {
                    flangebcd.setError("Fill this");
                    flangebcd.requestFocus();}
                if (flangeid.getText().toString().isEmpty()) {
                    flangeid.setError("Fill this");
                    flangeid.requestFocus();}
                if (flangeod.getText().toString().isEmpty()) {
                    flangeod.setError("Fill this");
                    flangeod.requestFocus();}
                if (gasketcontactarea.getText().toString().isEmpty()) {
                    gasketcontactarea.setError("Fill this");
                    gasketcontactarea.requestFocus();}
                if (gasketid.getText().toString().isEmpty()) {
                    gasketid.setError("Fill this");
                    gasketid.requestFocus();}
                if (gasketod.getText().toString().isEmpty()) {
                    gasketod.setError("Fill this");
                    gasketod.requestFocus();}
                if (temp.getText().toString().isEmpty()){
                    temp.setError("Fill this");
                    temp.requestFocus();}
                if (pressure.getText().toString().isEmpty()){         //Checking if the input is empty or not
                    pressure.setError("Fill this");                  //Shows error message in the writing area for filling up
                    pressure.requestFocus();}                        //returns cursor to the empty writing portal

                try {                                                       //tries this, if finds any error, goes to catch
                    calculateflange();                                      //Calling function createpdfwrapper
                } catch (FileNotFoundException e) {
                    e.printStackTrace();                                    //Helps in finding errors while in debugging mode
                } catch (DocumentException e) {
                    e.printStackTrace();                                    //Helps in finding errors while in debugging mode
                }
            }
        });
    }




    //Function #2---This is a function which consists of all calculation steps
    private void calculateflange() throws FileNotFoundException,DocumentException{
        //This is getting numerical values of inputs to the new variables
        String Pring = (pressure.getText().toString());
        String TEring = (temp.getText().toString());
        String FAring = (force.getText().toString());
        String MEring = (moment.getText().toString());
        String Godring = (gasketod.getText().toString());
        String Gidring = (gasketid.getText().toString());
        String Gcring = (gasketcontactarea.getText().toString());
        String Aring = (flangeod.getText().toString());
        String Bring = (flangeid.getText().toString());
        String Cring = (flangebcd.getText().toString());
        String tring = (flangethk.getText().toString());
        String goring = (hubsmall.getText().toString());
        String g1ring = (hublarge.getText().toString());
        String hring = (hublength.getText().toString());
        String KRring = (rigidity.getText().toString());
        String nbring = (boltnum.getText().toString());
        String dbring = (boltdia.getText().toString());


        //This is assigning variables which are used in calculations
        float N;
        float bo;
        float b;
        float G;
        float Cul=1;
        float Cus=1;
        float Wo;
        float Am;
        float Sbo=23000;
        float Sbg=23000;
        float Sfo=20000;
        float Sfg=20000;
        float Sno=17100;
        float Sng=17100;
        float m=3;
        float y=10000;
        float Ab;
        float Wg;
        float SF;


        float ho;
        float F;
        float V;
        float f;
        float hD;
        float AA;
        float BB;
        float CC;
        float DDG;
        float Mo;
        float Mg;
        float Fs=1;
        float alpha=1;
        float SHo;
        float SHg;
        float SRo;
        float SRg;
        float STo1;
        float STg1;
        float STo2=0;
        float STg2=0;
        int Eyo=28800000;
        int Eyg=29500000;



        int PT;     //0 for internal pressure     1 for reverse condition
        int GT;     //0 for NSEG gasket           1 for Self energized gasket
        int FT;     //0 for integral type         1 for loose type flange
        if(ptype.isChecked()){ PT = 1; }else{ PT=0; }
        if(gtype.isChecked()){ GT = 1; }else{ GT=0; }
        if(ftype.isChecked()){ FT = 1; }else{ FT=0; }



        if(!Pring.isEmpty() && !TEring.isEmpty() && !Godring.isEmpty() && !Gidring.isEmpty() && !Gcring.isEmpty() && !Aring.isEmpty() && !Bring.isEmpty() && !Cring.isEmpty() && !tring.isEmpty() && !goring.isEmpty() && !g1ring.isEmpty() && !hring.isEmpty() && !nbring.isEmpty() && !dbring.isEmpty()){
            Toast.makeText(getApplicationContext(),"Calculated !!!",Toast.LENGTH_LONG).show();

            int P = Integer.parseInt(Pring);
            int TE = Integer.parseInt(TEring);
            float FA = 0;
            if(!FAring.isEmpty()){FA = Float.parseFloat(FAring);}
            float ME = 0;
            if(!MEring.isEmpty()){ME = Float.parseFloat(MEring);}
            float God = Float.parseFloat(Godring);
            float Gid = Float.parseFloat(Gidring);
            float Gc = Float.parseFloat(Gcring);
            float A = Float.parseFloat(Aring);
            float B = Float.parseFloat(Bring);
            float C = Float.parseFloat(Cring);
            float t = Float.parseFloat(tring);
            float go = Float.parseFloat(goring);
            float g1 = Float.parseFloat(g1ring);
            float h = Float.parseFloat(hring);
            float KR;
            if(KRring.isEmpty()){ if(FT==0){ KR = (float)0.3; }else{ KR = (float)0.2; } }else{ KR = Float.parseFloat(KRring); }
            int nb = Integer.parseInt(nbring);
            float db = Float.parseFloat(dbring);


            DecimalFormat dff =new DecimalFormat("#.####");
            //dff.roundingMode = RoundingMode.CEILING;




            N = (God-Gid)/2;
            bo = N/2;

            if(bo <= 0.25){
                b = bo;
                G = (God+Gid)/2;}
            else{
                b = (float)(0.5*Cul*Math.sqrt(bo/Cul));
                G = Gc - 2 * b;}


            //if gasket is NSEG
            if(GT==0){
                Wo = (float)(0.785*G*G*P + 2*b*Math.PI*G*m*P);
                Am = (float)(Math.max(((Wo+FA+4*ME/G)/Sbo),(Math.PI*b*G*Cus*y/Sbg)));}
            else{     //if gasket is SEG
                Wo = (float)(0.785*G*G*P);
                Am = Math.max(((Wo+FA+4*ME/G)/Sbo),0);}

            Ab = (float)(nb*Math.PI*db*db/4);

            Wg = (Am+Ab)*Sbg/2;

            SF = Ab/Am;



            //flange calculation will start from here
            //***************************************
            float K = A/B;

            //if the flange is reverse
            if(PT == 1){ alpha = (float)((1/Math.pow(K,2)) * (1+0.668*(K+1)/ (1/(K-1)*(0.66845+5.7169*Math.pow(K,2)*Math.log10(K)/(Math.pow(K,2)-1)))));}

            float Y = (float)(alpha/(K-1)* (0.66845+5.7169*Math.pow(K,2)*Math.log10(K)/(Math.pow(K,2)-1)));

            float Z = (float)((Math.pow(K,2)+1)/(Math.pow(K,2)-1));

            float T = (float)((Math.pow(K,2)* (1+8.55246*Math.log10(K))-1) / ((1.0472+1.9448*Math.pow(K,2))*(K-1)));
            if(PT == 1) {
                T = (float)(alpha*T*(Z+0.3)/(Z-0.3));
            }

            float U = (float)(alpha*(Math.pow(K,2)* (1+8.55246*Math.log10(K))-1) / (1.36136*(Math.pow(K,2)-1)*(K-1)));



            if(PT == 0) {
                ho = (float)(Math.sqrt(B*go));
            }else{
                ho = (float)(Math.sqrt(A*go));
            }



            float Xg = g1/go;

            float Xh = h/ho;


            //for integral or reverse integral
            if(FT == 0){

                F = (float)(0.897697-0.297012*Math.log(Xg)+0.0095257*Math.log(Xh)+0.123586*Math.pow(Math.log(Xg),2)+0.035858*Math.pow(Math.log(Xh),2)-0.194422*Math.log(Xg)*Math.log(Xh)-0.0181259*Math.pow(Math.log(Xg),3)+0.012936*Math.pow(Math.log(Xh),3)-0.0377693*Math.log(Xg)*Math.pow(Math.log(Xh),2)+0.0273791*Math.log(Xh)*Math.pow(Math.log(Xg),2));


                if (0.1 <= Xh && Xh <= 0.5) {
                    V = (float)(0.500244 + 0.227914 / Xg - 1.87071 * Xh - 0.34441 / Math.pow(Xg, 2) + 2.49189 * Math.pow(Xh, 2) + 0.873446 * Xh / Xg + 0.189953 / Math.pow(Xg, 3) - 1.06082 * Math.pow(Xh, 3) - 1.4997 * Math.pow(Xh, 2) / Xg + 0.719413 * Xh / Math.pow(Xg, 2));
                } else {
                    V = (float)(0.0144868 - 0.135977 / Xg - 0.0461919 / Xh + 0.560718 / Math.pow(Xg, 2) + 0.0529829 / Math.pow(Xh, 2) + 0.244313 / (Xg * Xh) + 0.113929 / Math.pow(Xg, 3) - 0.00928265 / Math.pow(Xh, 3) - 0.0266293 / (Xg * Math.pow(Xh, 2)) - 0.217008 / (Math.pow(Xg, 2) * Xh));
                }


                f = (float)(Math.max(1,(0.0927779-0.0336633*Xg+0.964176*Math.pow(Xg,2)+0.0566286*Xh+0.347074*Math.pow(Xh,2)-4.18699*Math.pow(Xh,3))/(1-0.00596093*Xg+1.62904*Xh+3.49329*Math.pow(Xh,2)+1.39052*Math.pow(Xh,3))));

            }else{

                F = (float)((0.941074+0.176139*Math.log(Xg)-0.188556*Math.log(Xh)+0.0689847*Math.pow(Math.log(Xg),2)+0.523798*Math.pow(Math.log(Xh),2)-0.513894*Math.log(Xg)*Math.log(Xh))/(1+0.379392*Math.log(Xg)+0.18452*Math.log(Xh)-0.00605208*Math.pow(Math.log(Xg),2)-0.00358934*Math.pow(Math.log(Xh),2)+0.110179*Math.log(Xg)*Math.log(Xh)));


                if(0.1 <= Xh && Xh <= 0.25){
                    V = (float)(Math.exp(6.57683 - 0.115516 * Xg + 1.39499 * Math.sqrt(Xg) * Math.log(Xg) + 0.30734 * Math.pow(Math.log(Xg),2) - 8.30849 * Math.sqrt(Xg) + 2.62307 * Math.log(Xg) + 0.239498 * Xh * Math.log(Xh) - 2.96125 * Math.log(Xh) + 0.0007035052 / Xh));
                }else if(0.25 < Xh && Xh <= 0.5){
                    V = (float)(1.56323-1.80696*Math.log(Xg)-1.33458/Xh+0.276415*Math.pow(Math.log(Xg),2)+0.417135/Math.pow(Xh,2)+1.39511*Math.log(Xg)/Xh+0.0137129*Math.pow(Math.log(Xg),3)+0.0943597/Math.pow(Math.log(Xh),3)-0.402096*Math.log(Xg)/Math.pow(Xh,2)-0.101619*Math.pow(Math.log(Xg),2)/Xh);
                }else if( 0.5 < Xh && Xh <= 1){
                    V = (float)(-0.0213643-0.0763597/Xg+0.10299/Xh+0.725776/Math.pow(Xg,2)-0.160603/Math.pow(Xh,2)-0.0918061/(Xg*Xh)+0.472277/Math.pow(Xg,3)+0.087353/Math.pow(Xh,3)+0.527487/(Xg*Math.pow(Xh,2))-0.980209/(Math.pow(Xg,2)*Xh));
                }else{
                    V = (float)(0.00796687-0.220518/Xg+0.0602652/Xh+0.619818/Math.pow(Xg,2)-0.223212/Math.pow(Xh,2)+0.42192/(Xg*Xh)+0.0950195/Math.pow(Xg,3)+0.209813/Math.pow(Xh,3)-0.158821/(Xg*Math.pow(Xh,2))-0.242056/(Math.pow(Xg,2)*Xh));
                }

                f = 1;

            }



            float e = F/ho;

            float d = (float)(U*Math.pow(go,2)*ho/V);

            float L = (float)((t*e+1)/T + Math.pow(t,3)/d);

            float H = (float)(0.785*Math.pow(G,2)*P);

            float HD = (float)(0.785*Math.pow(B,2)*P);

            float HT = H-HD;

            float HG = Wo-H;


            //integral flange
            if(FT==0 && PT==0){
                hD = (C-B-g1)/2;
                //reverse integral
            }else if(FT==0 && PT==1){
                hD =  (C-B+g1-2*go)/2;
            }else{
                hD = (C-B)/2;
            }

            float hG = (C-G)/2;

            float hT = (float)(0.5*(C-(B+G)/2));

            float Bs = C/nb;

            float Bsc = (float)(Math.max(1,Math.sqrt(Bs/(2*db+t))));

            float I = (float)(0.0874*L*Math.pow(go,2)*ho*B/V);

            float AR = (float)(0.5*(A-B));

            float Gavg = (float)(0.5*(go+g1));

            if(t >= Gavg){
                AA = AR;
                BB = t;
                CC = h;
                DDG = Gavg;
            }else{
                AA = h+t;
                BB = Gavg;
                CC = AR-Gavg;
                DDG = t;
            }

            float KAB = (float)((AA*Math.pow(BB,3)) * (0.33333-0.21*BB/AA*(1-0.08333*Math.pow(BB/AA,4))));

            float KCD = (float)((CC * Math.pow(DDG, 3)) * (0.33333 - (0.105 * (DDG / CC) * (1 -0.005208* Math.pow((DDG / CC), 4)))));

            float Ip = KAB+KCD;

            float Moe = (float)(4*ME*(I/(0.3846*Ip+1))*(hD/(C-2*hD))+FA*hD);


            //internal pressure
            if(PT==0){
                Mo = (float)(Math.abs(((HD*hD+HT*hT+HG*hG)*Bsc+Moe)*Fs));
                Mg = Wg*(C-G)*Bsc*Fs/2;
            }else{
                Mo = (float)(Math.abs((HD*(hD-hG)+HT*(hT-hG)+Moe)*Fs));
                Mg = Wg*hG*Fs;
            }


            //Flange Stress Equations
            //***********************
            SHo = (float)(f*Mo/(L*Math.pow(g1,2)*B));
            SHg = (float)(f*Mg/(L*Math.pow(g1,2)*B));

            SRo = (float)((1.33*t*e+1)*Mo/(L*Math.pow(t,2)*B));
            SRg = (float)((1.33*t*e+1)*Mg/(L*Math.pow(t,2)*B));


            //integral or loose (internal pressure)
            if(PT==0){
                STo1 = (float)(Y*Mo/(Math.pow(t,2)*B)-Z*SRo);
                STg1 = (float)(Y*Mg/(Math.pow(t,2)*B)-Z*SRg);
            }else{
                STo1 = (float)(Y*Mo/(Math.pow(t,2)*B)-Z*SRo *(0.67*t*e+1)/(1.33*t*e+1));
                STg1 = (float)(Y*Mg/(Math.pow(t,2)*B)-Z*SRg *(0.67*t*e+1)/(1.33*t*e+1));

                STo2 = (float)(Mo/(Math.pow(t,2)*B) * (Y -2*Math.pow(K,2)*(0.67*t*e+1)/((Math.pow(K,2)-1)*L)));
                STg2 = (float)(Mg/(Math.pow(t,2)*B) * (Y -2*Math.pow(K,2)*(0.67*t*e+1)/((Math.pow(K,2)-1)*L)));
            }


            //Flange Rigidity Criterion
            //*************************
            float Jo = (float)(52.14*V*Mo/(L*Eyo*Math.pow(go,2)*KR*ho));
            float Jg = (float)(52.14*V*Mg/(L*Eyg*Math.pow(go,2)*KR*ho));



            //resultN.setText("N: " + dff.format(N) );


            createPdfWrapper(P,TE,FA,ME,God,Gid,Gc,A,B,C,t,go,g1,h,KR,nb,db,N,bo,b,G,Wo,Am,Ab,Wg,SF,PT,GT,FT,m,y,Sbo,Sbg,Sfo,Sfg,Sno,Sng,
                    ho,F,V,hD,AA,BB,CC,DDG,Mo,Mg,alpha,SHo,SHg,SRo,SRg,STo1,STg1,STo2,STg2,Eyo,Eyg,
                    K,Z,Y,T,U,Xg,Xh,e,d,L,f,H,HD,HT,HG,hG,hT,Bs,Bsc,I,AR,Gavg,KAB,KCD,Ip,Moe,Jo,Jg);

        }else{
            Toast.makeText(getApplicationContext(),"Please, fill all the boxes!",Toast.LENGTH_LONG).show();
        }

    }



    //This is function #3-----This function checks for permission to write, if permission is granted it writes, if not, asks for permission
    private void createPdfWrapper(int P5, int TE5, float FA5, float ME5, float God5, float Gid5, float Gc5, float A5, float B5, float C5, float t5, float go5, float g15, float h5, float KR5, int nb5, float db5, float N5, float bo5, float b5, float G5, float Wo5, float Am5, float Ab5, float Wg5, float SF5, int PT5, int GT5, int FT5, float m5, float y5, float Sbo5, float Sbg5, float Sfo5, float Sfg5, float Sno5, float Sng5,
                                  float ho5, float F5, float V5, float hD5, float AA5, float BB5, float CC5, float DDG5, float Mo5, float Mg5, float alpha5, float SHo5, float SHg5, float SRo5, float SRg5, float STo15, float STg15, float STo25, float STg25, float Eyo5, float Eyg5,
                                  float K5, float Z5, float Y5, float T5, float U5, float Xg5, float Xh5, float e5, float d5, float L5, float f5, float H5, float HD5, float HT5, float HG5, float hG5, float hT5, float Bs5, float Bsc5, float I5, float AR5, float Gavg5, float KAB5, float KCD5, float Ip5, float Moe5, float Jo5, float Jg5) throws FileNotFoundException,DocumentException{
        int pertowrite = ActivityCompat.checkSelfPermission(this, Manifest.permission.WRITE_EXTERNAL_STORAGE);
        if (pertowrite != PackageManager.PERMISSION_GRANTED) {
            if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.M) {
                requestPermissions(new String[]{Manifest.permission.WRITE_EXTERNAL_STORAGE},REQUEST_CODE_ASK_PERMISSIONS); }
        }else{
            createPdf(P5,TE5,FA5,ME5,God5,Gid5,Gc5,A5,B5,C5,t5,go5,g15,h5,KR5,nb5,db5,N5,bo5,b5,G5,Wo5,Am5,Ab5,Wg5,SF5,PT5,GT5,FT5,m5,y5,Sbo5,Sbg5,Sfo5,Sfg5,Sno5,Sng5,
                    ho5,F5,V5,hD5,AA5,BB5,CC5,DDG5,Mo5,Mg5,alpha5,SHo5,SHg5,SRo5,SRg5,STo15,STg15,STo25,STg25,Eyo5,Eyg5,
                    K5,Z5,Y5,T5,U5,Xg5,Xh5,e5,d5,L5,f5,H5,HD5,HT5,HG5,hG5,hT5,Bs5,Bsc5,I5,AR5,Gavg5,KAB5,KCD5,Ip5,Moe5,Jo5,Jg5);
        }
    }




    //This is function #4-----This function asks yes no for any confirmation
    private void showMessageOKCancel(String message, DialogInterface.OnClickListener okListener) {
        new AlertDialog.Builder(this)
                .setMessage(message)
                .setPositiveButton("Yup", okListener)
                .setNegativeButton("Nope", null)
                .create()
                .show();
    }



    //This is function #5-----This function makes contents for PDF file, writes and make a pdf file in instructed directory
    private void createPdf(int P6, int TE6, float FA6, float ME6, float God6, float Gid6, float Gc6, float A6, float B6, float C6, float t6, float go6, float g16, float h6, float KR6, int nb6, float db6, float N6, float bo6, float b6, float G6, float Wo6, float Am6, float Ab6, float Wg6, float SF6, int PT6, int GT6, int FT6, float m6, float y6, float Sbo6, float Sbg6, float Sfo6, float Sfg6, float Sno6, float Sng6,
                           float ho6, float F6, float V6, float hD6, float AA6, float BB6, float CC6, float DDG6, float Mo6, float Mg6, float alpha6, float SHo6, float SHg6, float SRo6, float SRg6, float STo16, float STg16, float STo26, float STg26, float Eyo6, float Eyg6,
                           float K6, float Z6, float Y6, float T6, float U6, float Xg6, float Xh6, float e6, float d6, float L6, float f6, float H6, float HD6, float HT6, float HG6, float hG6, float hT6, float Bs6, float Bsc6, float I6, float AR6, float Gavg6, float KAB6, float KCD6, float Ip6, float Moe6, float Jo6, float Jg6) throws FileNotFoundException, DocumentException {         //function that is to be called for creating pdf
        File docsFolder = new File(Environment.getExternalStorageDirectory() + "/FlangeCALC");      //variable assigned as docsfolder representing the file directory
        if (!docsFolder.exists()) {
            docsFolder.mkdir();          //instruction to make-directory as this docsfolder variable is assigned to
        }

        String timeStamp = new SimpleDateFormat("MM_dd_yyyy___HHmmss").format(new Date());
        pdfFile = new File(docsFolder.getAbsolutePath(),"Flange_Calculator_"+ timeStamp + ".pdf");    //requesting to make new pdf file, try to rename file with time dependent
        OutputStream output = new FileOutputStream(pdfFile);
        Document document = new Document();
        PdfWriter.getInstance(document, output);



        DecimalFormat dff =new DecimalFormat("#.####");
        float FA7 = Float.parseFloat(dff.format(FA6));
        float ME7 = Float.parseFloat(dff.format(ME6));
        //float xyz = Float.parseFloat(dff.format(55555));
        float God7 = Float.parseFloat(dff.format(God6));
        float Gid7 = Float.parseFloat(dff.format(Gid6));
        float Gc7 = Float.parseFloat(dff.format(Gc6));
        //float xyz = Float.parseFloat(dff.format(55555));
        //float xyz = Float.parseFloat(dff.format(55555));
        //float xyz = Float.parseFloat(dff.format(55555));
        //float xyz = Float.parseFloat(dff.format(55555));
        //float xyz = Float.parseFloat(dff.format(55555));
        float A7 = Float.parseFloat(dff.format(A6));
        float B7 = Float.parseFloat(dff.format(B6));
        float C7 = Float.parseFloat(dff.format(C6));
        float t7 = Float.parseFloat(dff.format(t6));
        float go7 = Float.parseFloat(dff.format(go6));
        float g17 = Float.parseFloat(dff.format(g16));
        float h7 = Float.parseFloat(dff.format(h6));
        float KR7 = Float.parseFloat(dff.format(KR6));
        float db7 = Float.parseFloat(dff.format(db6));
        //float xyz = Float.parseFloat(dff.format(xyzz));



        float N7 = Float.parseFloat(dff.format(N6));
        float bo7 = Float.parseFloat(dff.format(bo6));
        float b7 = Float.parseFloat(dff.format(b6));
        float G7 = Float.parseFloat(dff.format(G6));
        float Wo7 = Float.parseFloat(dff.format(Wo6));
        float Am7 = Float.parseFloat(dff.format(Am6));
        float Ab7 = Float.parseFloat(dff.format(Ab6));
        float Wg7 = Float.parseFloat(dff.format(Wg6));
        float SF7 = Float.parseFloat(dff.format(SF6));
        float PT7 = Float.parseFloat(dff.format(PT6));
        float GT7 = Float.parseFloat(dff.format(GT6));
        float FT7 = Float.parseFloat(dff.format(FT6));
        float m7 = Float.parseFloat(dff.format(m6));
        float y7 = Float.parseFloat(dff.format(y6));
        float ho7 = Float.parseFloat(dff.format(ho6));
        float F7 = Float.parseFloat(dff.format(F6));
        float V7 = Float.parseFloat(dff.format(V6));
        float hD7 = Float.parseFloat(dff.format(hD6));
        float AA7 = Float.parseFloat(dff.format(AA6));
        float BB7 = Float.parseFloat(dff.format(BB6));
        float CC7 = Float.parseFloat(dff.format(CC6));
        float DDG7 = Float.parseFloat(dff.format(DDG6));
        float Mo7 = Float.parseFloat(dff.format(Mo6));
        float Mg7 = Float.parseFloat(dff.format(Mg6));
        float alpha7 = Float.parseFloat(dff.format(alpha6));
        float SHo7 = Float.parseFloat(dff.format(SHo6));
        float SHg7 = Float.parseFloat(dff.format(SHg6));
        float SRo7 = Float.parseFloat(dff.format(SRo6));
        float SRg7 = Float.parseFloat(dff.format(SRg6));
        float STo17 = Float.parseFloat(dff.format(STo16));
        float STg17 = Float.parseFloat(dff.format(STg16));
        float STo27 = Float.parseFloat(dff.format(STo26));
        float STg27 = Float.parseFloat(dff.format(STg26));
        float Eyo7 = Float.parseFloat(dff.format(Eyo6));
        float Eyg7 = Float.parseFloat(dff.format(Eyg6));
        float K7 = Float.parseFloat(dff.format(K6));
        float Z7 = Float.parseFloat(dff.format(Z6));
        float Y7 = Float.parseFloat(dff.format(Y6));
        float T7 = Float.parseFloat(dff.format(T6));
        float U7 = Float.parseFloat(dff.format(U6));
        float Xg7 = Float.parseFloat(dff.format(Xg6));
        float Xh7 = Float.parseFloat(dff.format(Xh6));
        float e7 = Float.parseFloat(dff.format(e6));
        float d7 = Float.parseFloat(dff.format(d6));
        float L7 = Float.parseFloat(dff.format(L6));
        float f7 = Float.parseFloat(dff.format(f6));
        float H7 = Float.parseFloat(dff.format(H6));
        float HD7 = Float.parseFloat(dff.format(HD6));
        float HT7 = Float.parseFloat(dff.format(HT6));
        float HG7 = Float.parseFloat(dff.format(HG6));
        float hG7 = Float.parseFloat(dff.format(hG6));
        float hT7 = Float.parseFloat(dff.format(hT6));
        float Bs7 = Float.parseFloat(dff.format(Bs6));
        float Bsc7 = Float.parseFloat(dff.format(Bsc6));
        float I7 = Float.parseFloat(dff.format(I6));
        float AR7 = Float.parseFloat(dff.format(AR6));
        float Gavg7 = Float.parseFloat(dff.format(Gavg6));
        float KAB7 = Float.parseFloat(dff.format(KAB6));
        float KCD7 = Float.parseFloat(dff.format(KCD6));
        float Ip7 = Float.parseFloat(dff.format(Ip6));
        float Moe7 = Float.parseFloat(dff.format(Moe6));
        float Jo7 = Float.parseFloat(dff.format(Jo6));
        float Jg7 = Float.parseFloat(dff.format(Jg6));



        int areacheck;
        if(Am7 <= Ab7){areacheck = 1;}else{areacheck = 2;}


        int Jocheck;
        if(Jo7 <= 1){Jocheck = 1;}else{Jocheck = 2;}


        int Jgcheck;
        if(Jg7 <= 1){Jgcheck = 1;}else{Jgcheck = 2;}


        int stressOC = 1;
        if(PT7 == 0){               ///////For integral and loose type flange
            //First condition
            if(FT7 == 0){           ///////for integral type
                if(SHo7 > Math.min(1.5*Sfo6,2.5*Sno6)){ stressOC = 2; }
            }else{                  ///////for loose type
                if(SHo7 > 1.5*Sfo6){ stressOC = 2; } }
            //Second condition
            if(SRo7 > Sfo6){ stressOC = 2; }
            //Third condition
            if(STo17 > Sfo6){ stressOC = 2; }
            //Fourth condition
            if((SHo7+SRo7)/2 > Sfo6){ stressOC = 2; }
            //Fifth condition
            if((SHo7+STo17)/2 > Sfo6){ stressOC = 2; }
        }else{                      ///////For reverse integral and loose type flange
            //First condition
            if(SHo7 > 1.5*Sfo6){ stressOC = 2; }
            //Second condition
            if(SRo7 > Sfo6){ stressOC = 2; }
            //Third condition
            if(STo17 > Sfo6){ stressOC = 2; }
            //Fourth condition
            if((SHo7+SRo7)/2 > Sfo6){ stressOC = 2; }
            //Fifth condition
            if((SHo7+STo17)/2 > Sfo6){ stressOC = 2; }
            //Sixth condition
            if(STo27 > Sfo6){ stressOC = 2; }
        }


        int stressGSC = 1;
        if(PT7 == 0){               ///////For integral and loose type flange
            //First condition
            if(FT7 == 0){           ///////for integral type
                if(SHg7 > Math.min(1.5*Sfg6,2.5*Sng6)){ stressGSC = 2; }
            }else{                  ///////for loose type
                if(SHg7 > 1.5*Sfg6){ stressGSC = 2; } }
            //Second condition
            if(SRg7 > Sfg6){ stressGSC = 2; }
            //Third condition
            if(STg17 > Sfg6){ stressGSC = 2; }
            //Fourth condition
            if((SHg7+SRg7)/2 > Sfg6){ stressGSC = 2; }
            //Fifth condition
            if((SHg7+STg17)/2 > Sfg6){ stressGSC = 2; }
        }else{                      ///////For reverse integral and loose type flange
            //First condition
            if(SHg7 > 1.5*Sfg6){ stressGSC = 2; }
            //Second condition
            if(SRg7 > Sfg6){ stressGSC = 2; }
            //Third condition
            if(STg17 > Sfg6){ stressGSC = 2; }
            //Fourth condition
            if((SHg7+SRg7)/2 > Sfg6){ stressGSC = 2; }
            //Fifth condition
            if((SHg7+STg17)/2 > Sfg6){ stressGSC = 2; }
            //Sixth condition
            if(STg27 > Sfg6){ stressGSC = 2; }
        }







        String parameterinput[] = {"Design Pressure","Design Temperature","External Axial Force","External Bending Moment","Gasket Material","Gasket Outer Diameter","Gasket Inside Diameter","Gasket Contact Area Diameter","Gasket Type Selected","Flange Material","Pipe Material","Bolt Material","Flange Type Selected","Flange Outer Diameter","Flange Inside Diameter","Flange Bolt Circle Diameter","Flange Thickness","Small Hub Thickness","Large Hub Thickness","Flange Hub Length","Flange Rigidity Factor","Number of Bolts in Flange","Bolt Diameter"};
        String notationinput[] = {"P","T","FA","ME","G-Material","God","Gid","Gc","G-Type","F-Material","P-Material","B-Material","F-Type","A","B","C","t","go","g1","h","KR","nb","db"};
        float valueinput[] = {P6,TE6,FA7,ME7,55555,God7,Gid7,Gc7,55555,55555,55555,55555,55555,A7,B7,C7,t7,go7,g17,h7,KR7,nb6,db7};

        PdfPTable ptinput = new PdfPTable(new float[] {5,50,20,25});
        ptinput.getDefaultCell().setHorizontalAlignment(Element.ALIGN_CENTER);
        ptinput.addCell("SN");
        ptinput.addCell("Input Parameters");
        ptinput.addCell("Notations");
        ptinput.addCell("Values");
        ptinput.setHeaderRows(1);
        PdfPCell[] cellinput = ptinput.getRow(0).getCells();
        for (int j=0; j<cellinput.length; j++){
            cellinput[j].setBackgroundColor(BaseColor.LIGHT_GRAY);
        }
        for (int i=1; i<24; i++){
            ptinput.addCell(""+i);
            ptinput.addCell(""+parameterinput[i-1]);
            ptinput.addCell(""+notationinput[i-1]);
            ptinput.addCell(""+valueinput[i-1]);
        }




        String parameterresult[] = {"Effective Gasket Contact Width","Gasket Load Reaction Diameter","Minimum Required Bolt Area","Design Bolt Load (OC)","Design Bolt Load (GSC)","Bolt Safety Factor","Total Hydrostatic End Force","Bending Moment of Inertia FCS","Polar Moment of Inertia FCS","Flange Design Moment (OC)","Flange Design Moment (GSC)","Flange Stress Acceptance (OC)","Flange Stress Acceptance (GSC)","Flange Rigidity (OC)","Flange Rigidity (GSC)"};
        String notationresult[] = {"b","G","Am","Wo","Wg","SF","H","I","Ip","Mo","Mg","N/A","N/A","Jo","Jg"};
        float valueresult[] = {b7,G7,Am7,Wo7,Wg7,SF7,H7,I7,Ip7,Mo7,Mg7,-1,-1,Jo7,Jg7};
        int passcondition[] = {0,0,areacheck,0,0,areacheck,0,0,0,0,0,stressOC,stressGSC,Jocheck,Jgcheck};

        PdfPTable ptresult = new PdfPTable(new float[] {5,44,14,20,17});
        ptresult.getDefaultCell().setHorizontalAlignment(Element.ALIGN_CENTER);
        ptresult.addCell("SN");
        ptresult.addCell("Result Parameters");
        ptresult.addCell("Notations");
        ptresult.addCell("Values");
        ptresult.addCell("Remarks");
        ptresult.setHeaderRows(1);
        PdfPCell[] cellresult = ptresult.getRow(0).getCells();
        for (int j=0; j<cellresult.length; j++){
            cellresult[j].setBackgroundColor(BaseColor.LIGHT_GRAY); }
        for (int i=1; i<16; i++){
            ptresult.addCell(""+i);
            ptresult.addCell(""+parameterresult[i-1]);
            ptresult.addCell(""+notationresult[i-1]);
            if(valueresult[i-1]>0){
                ptresult.addCell(""+valueresult[i-1]);
            }else{
                ptresult.addCell("N/A"); }

            if(passcondition[i-1] == 0){
                ptresult.addCell("N/A");
            }else if(passcondition[i-1] == 1){
                ptresult.addCell("Passed");
            }else{
                ptresult.addCell("Not Passed"); } }





        document.open();       //opens file and start writing in it
        document.add(new Paragraph("Flange Calculator - 2019  ©Sushil Champ"));
        document.add(new Paragraph("Report created on: "+new Date()));
        //document.add(new Paragraph("****************************************************************************************************************\n"));
        document.add(new Paragraph("****************************************************************************************************************\n"));
        document.add(new Paragraph("1.   Inputs Provided:\n\n"));
        document.add(new Paragraph("Following are the inputs provided:\n"));
        document.add(new Paragraph("*All calculations are done in US Customary Units.\n\n"));
        document.add(ptinput);


        int i = document.getPageNumber();


        document.newPage();
        document.add(new Paragraph("Flange Calculator - 2019  ©Sushil Champ"));
        document.add(new Paragraph("Report created on: "+new Date()));
        document.add(new Paragraph("****************************************************************************************************************\n"));
        document.add(new Paragraph("2.   Calculation and Results:\n"));
        document.add(new Paragraph("Design calculations for bolt loads and flange loads are performed as follow:\n\n"));
        document.add(new Paragraph("A)   Design for Bolt Loads:\n\n"));

        document.add(new Paragraph("Gasket width: N\n"));
        document.add(new Paragraph("= (God-Gid)/2\n"));
        document.add(new Paragraph("= (" +God7+ "-" +Gid7+ ")/2\n"));
        document.add(new Paragraph("= " +N7+ "  in\n\n"));


        document.add(new Paragraph("Basic gasket seating width: bo\n"));
        document.add(new Paragraph("= N/2\n"));
        document.add(new Paragraph("= " +N7+ "/2\n"));
        document.add(new Paragraph("= " +bo7+ "  in\n\n"));


        document.add(new Paragraph("Effective gasket seating width: b\n"));
        if(bo6 <= 0.25){
            document.add(new Paragraph("= bo\n"));
        }else{
            document.add(new Paragraph("= 0.5*Cul*sqrt(bo/Cul)\n"));
            document.add(new Paragraph("= 0.5*1*sqrt(" +bo7+ "/1)\n")); }
        document.add(new Paragraph("= " +b7+ "  in\n\n"));


        document.add(new Paragraph("Location of gasket reaction: G\n"));
        if(bo6 <= 0.25){
            document.add(new Paragraph("= (God+Gid)/2\n"));
            document.add(new Paragraph("= (" +God7+ "+" +Gid7+ ")/2\n"));
        }else{
            document.add(new Paragraph("= Gc-2*b\n"));
            document.add(new Paragraph("= " +Gc7+ "-2*" +b7+ "\n")); }
        document.add(new Paragraph("= " +G7+ "  in\n\n"));


        document.add(new Paragraph("Design bolt load for the operating condition: Wo\n"));
        if(GT6 == 0){
            document.add(new Paragraph("= 0.785*G^2*P + 2*b*PI*G*m*P\n"));
            document.add(new Paragraph("= 0.785*" +G7+ "^2*" +P6+ " + 2*" +b7+ "*PI*" +G7+ "*" +m7+ "*" +P6+ "\n"));
        }else{
            document.add(new Paragraph("= 0.785*G^2*P\n"));
            document.add(new Paragraph("= 0.785*" +G7+ "^2*" +P6+ "\n")); }
        document.add(new Paragraph("= " +Wo7+ "  lb.\n\n"));


        document.add(new Paragraph("Total minimum required cross-sectional area of bolts: Am\n"));
        if(GT6 == 0){
            document.add(new Paragraph("= max(((Wo+FA+4*ME/G)/Sbo),(PI*b*G*Cus*y/Sbg))\n"));
            document.add(new Paragraph("= max(((" +Wo7+ "+" +FA7+ "+4*" +ME7+ "/" +G7+ ")/" +Sbo6+ ") , (PI*" +b7+ "*" +G7+ "*1*" +y7+ "/" +Sbg6+ "))\n"));
        }else{
            document.add(new Paragraph("= max(((Wo+FA+4*ME/G)/Sbo),0)\n"));
            document.add(new Paragraph("= max(((" +Wo7+ "+" +FA7+ "+4*" +ME7+ "/" +G7+ ")/" +Sbo6+ "),0)\n")); }
        document.add(new Paragraph("= " +Am7+ "  in^2\n\n"));


        document.add(new Paragraph("Total actual cross-sectional area of bolts: Ab\n"));
        document.add(new Paragraph("= nb*PI*db^2/4\n"));
        document.add(new Paragraph("= " +nb6+ "*PI*" +db7+ "^2/4\n"));
        document.add(new Paragraph("= " +Ab7+ "  in^2\n\n"));


        document.add(new Paragraph("Design bolt load for the gasket seating condition: Wg\n"));
        document.add(new Paragraph("= (Am+Ab)*Sbg/2\n"));
        document.add(new Paragraph("= (" +Am7+ "+" +Ab7+ ")*" +Sbg6+ "/2\n"));
        document.add(new Paragraph("= " +Wg7+ "  lb.\n\n"));


        document.add(new Paragraph("Bolt safety factor: SF\n"));
        document.add(new Paragraph("= Ab/Am\n"));
        document.add(new Paragraph("= " +Ab7+ "/" +Am7+ "\n"));
        document.add(new Paragraph("= " +SF7+ "  \n\n"));
        if(SF7>=1){
            document.add(new Paragraph("Since" +SF7+ " >= 1, the bolt is safe under provided design conditions.\n"));
        }else{
            document.add(new Paragraph("Since" +SF7+ " < 1, the bolt is not safe under provided design conditions.\n"));
            document.add(new Paragraph("Increase bolt size & bolt number or reduce bolt loads.\n"));
        }


        document.newPage();
        document.add(new Paragraph("B)   Design for Flange Loads:\n\n"));
        document.add(new Paragraph("Flange Stress factors:\n\n"));


        document.add(new Paragraph("K     = A/B\n"));
        document.add(new Paragraph("= " +A7+ "/" +B7+ "\n"));
        document.add(new Paragraph("= " +K7+ "  \n\n"));


        if(PT6 == 1){
            document.add(new Paragraph("alpha     = (1/pow(K,2))*(1+0.668*(K+1)/(1/(K-1)*(0.66845+5.7169*pow(K,2)*log10(K)/(pow(K,2)-1))))\n"));
            document.add(new Paragraph("= (1/pow(" +K7+ ",2))*(1+0.668*(" +K7+ "+1)/(1/(" +K7+ "-1)*(0.66845+5.7169*pow(" +K7+ ",2)*log10(" +K7+ ")/(pow(" +K7+ ",2)-1))))\n"));
            document.add(new Paragraph("= " +alpha7+ "  \n\n")); }


        document.add(new Paragraph("Z     = (pow(K,2)+1)/(pow(K,2)-1)\n"));
        document.add(new Paragraph("= (pow(" +K7+ ",2)+1)/(pow(" +K7+ ",2)-1)\n"));
        document.add(new Paragraph("= " +Z7+ "  \n\n"));


        document.add(new Paragraph("Y     = alpha/(K-1)*(0.66845+5.7169*pow(K,2)*log10(K)/(pow(K,2)-1))\n"));
        document.add(new Paragraph("= " +alpha7+ "/(" +K7+ "-1)*(0.66845+5.7169*pow(" +K7+ ",2)*log10(" +K7+ ")/(pow(" +K7+ ",2)-1))\n"));
        document.add(new Paragraph("= " +Y7+ "  \n\n"));


        if(PT6 == 1){
            document.add(new Paragraph("T     = alpha*((Z+0.3)/(Z-0.3))*(pow(K,2)*(1+8.55246*log10(K))-1)/((1.0472+1.9448*pow(K,2))*(K-1))\n"));
            document.add(new Paragraph("= " +alpha7+ "*((" +Z7+ "+0.3)/(" +Z7+ "-0.3))*(pow(" +K7+ ",2)*(1+8.55246*log10(" +K7+ "))-1)/((1.0472+1.9448*pow(" +K7+ ",2))*(" +K7+ "-1))\n"));
        }else{
            document.add(new Paragraph("T     = alpha*(pow(K,2)*(1+8.55246*log10(K))-1)/((1.0472+1.9448*pow(K,2))*(K-1))\n"));
            document.add(new Paragraph("= " +alpha7+ "*(pow(" +K7+ ",2)*(1+8.55246*log10(" +K7+ "))-1)/((1.0472+1.9448*pow(" +K7+ ",2))*(" +K7+ "-1))\n"));}
        document.add(new Paragraph("= " +T7+ "  \n\n"));


        document.add(new Paragraph("U     = alpha*(pow(K,2)*(1+8.55246*log10(K))-1)/(1.36136*(pow(K,2)-1)*(K-1))\n"));
        document.add(new Paragraph("= " +alpha7+ "*(pow(" +K7+ ",2)*(1+8.55246*log10(" +K7+ "))-1)/(1.36136*(pow(" +K7+ ",2)-1)*(" +K7+ "-1))\n"));
        document.add(new Paragraph("= " +U7+ "  \n\n"));


        if(PT7 == 0){
            document.add(new Paragraph("ho     = sqrt(B*go)\n"));
            document.add(new Paragraph("= sqrt(" +B7+ "*" +go7+ ")\n"));
        }else{
            document.add(new Paragraph("ho     = sqrt(A*go)\n"));
            document.add(new Paragraph("= sqrt(" +A7+ "*" +go7+ ")\n"));}
        document.add(new Paragraph("= " +ho7+ "  \n\n"));


        document.add(new Paragraph("Xg     = g1/go\n"));
        document.add(new Paragraph("= " +g17+ "/" +go7+ "\n"));
        document.add(new Paragraph("= " +Xg7+ "  \n\n"));


        document.add(new Paragraph("Xh     = h/ho\n"));
        document.add(new Paragraph("= " +h7+ "/" +ho7+ "\n"));
        document.add(new Paragraph("= " +Xg7+ "  \n\n"));


        if(FT7 == 0){
            document.add(new Paragraph("F     = 0.897697-0.297012*log(Xg)+0.0095257*log(Xh)+0.123586*pow(log(Xg),2)+0.035858*pow(log(Xh),2)-0.194422*log(Xg)*log(Xh)-0.0181259*pow(log(Xg),3)+0.012936*pow(log(Xh),3)-0.0377693*log(Xg)*pow(log(Xh),2)+0.0273791*log(Xh)*pow(log(Xg),2)\n"));
            document.add(new Paragraph("= 0.897697-0.297012*log(" +Xg7+ ")+0.0095257*log(" +Xh7+ ")+0.123586*pow(log(" +Xg7+ "),2)+0.035858*pow(log(" +Xh7+ "),2)-0.194422*log(" +Xg7+ ")*log(" +Xh7+ ")-0.0181259*pow(log(" +Xg7+ "),3)+0.012936*pow(log(" +Xh7+ "),3)-0.0377693*log(" +Xg7+ ")*pow(log(" +Xh7+ "),2)+0.0273791*log(" +Xh7+ ")*pow(log(" +Xg7+ "),2)\n"));
        }else{
            document.add(new Paragraph("F     = (0.941074+0.176139*log(Xg)-0.188556*log(Xh)+0.0689847*pow(log(Xg),2)+0.523798*pow(log(Xh),2)-0.513894*log(Xg)*log(Xh))/(1+0.379392*log(Xg)+0.18452*log(Xh)-0.00605208*pow(log(Xg),2)-0.00358934*pow(log(Xh),2)+0.110179*log(Xg)*log(Xh))\n"));
            document.add(new Paragraph("= (0.941074+0.176139*log(" +Xg7+ ")-0.188556*log(" +Xh7+ ")+0.0689847*pow(log(" +Xg7+ "),2)+0.523798*pow(log(" +Xh7+ "),2)-0.513894*log(" +Xg7+ ")*log(" +Xh7+ "))/(1+0.379392*log(" +Xg7+ ")+0.18452*log(" +Xh7+ ")-0.00605208*pow(log(" +Xg7+ "),2)-0.00358934*pow(log(" +Xh7+ "),2)+0.110179*log(" +Xg7+ ")*log(" +Xh7+ "))\n")); }
        document.add(new Paragraph("= " +F7+ "  \n\n"));


        if(FT7 == 0){
            if(0.1 <= Xh7 && Xh7 <= 0.5){
                document.add(new Paragraph("V     = 0.500244+0.227914/Xg-1.87071*Xh-0.34441/pow(Xg,2)+2.49189*pow(Xh,2)+0.873446*Xh/Xg+0.189953/pow(Xg,3)-1.06082*pow(Xh,3)-1.4997*pow(Xh,2)/Xg+0.719413*Xh/pow(Xg,2)\n"));
                document.add(new Paragraph("= 0.500244+0.227914/" +Xg7+ "-1.87071*" +Xh7+ "-0.34441/pow(" +Xg7+ ",2)+2.49189*pow(" +Xh7+ ",2)+0.873446*" +Xh7+ "/" +Xg7+ "+0.189953/pow(" +Xg7+ ",3)-1.06082*pow(" +Xh7+ ",3)-1.4997*pow(" +Xh7+ ",2)/" +Xg7+ "+0.719413*" +Xh7+ "/pow(" +Xg7+ ",2)\n"));
            }else{
                document.add(new Paragraph("V     = 0.0144868-0.135977/Xg-0.0461919/Xh+0.560718/pow(Xg,2)+0.0529829/pow(Xh,2)+0.244313/(Xg*Xh)+0.113929/pow(Xg,3)-0.00928265/pow(Xh,3)-0.0266293/(Xg*pow(Xh,2))-0.217008/(pow(Xg,2)*Xh)\n"));
                document.add(new Paragraph("= 0.0144868-0.135977/" +Xg7+ "-0.0461919/" +Xh7+ "+0.560718/pow(" +Xg7+ ",2)+0.0529829/pow(" +Xh7+ ",2)+0.244313/(" +Xg7+ "*" +Xh7+ ")+0.113929/pow(" +Xg7+ ",3)-0.00928265/pow(" +Xh7+ ",3)-0.0266293/(" +Xg7+ "*pow(" +Xh7+ ",2))-0.217008/(pow(" +Xg7+ ",2)*" +Xh7+ ")\n")); }
        }else{
            if(0.1 <= Xh7 && Xh7 <= 0.25){
                document.add(new Paragraph("V     = exp(6.57683-0.115516*Xg+1.39499*sqrt(Xg)*log(Xg)+0.30734*pow(log(Xg),2)-8.30849*sqrt(Xg)+2.62307*log(Xg)+0.239498*Xh*log(Xh)-2.96125*log(Xh)+0.0007035052/Xh)\n"));
                document.add(new Paragraph("= exp(6.57683-0.115516*" +Xg7+ "+1.39499*sqrt(" +Xg7+ ")*log(" +Xg7+ ")+0.30734*pow(log(" +Xg7+ "),2)-8.30849*sqrt(" +Xg7+ ")+2.62307*log(" +Xg7+ ")+0.239498*" +Xh7+ "*log(" +Xh7+ ")-2.96125*log(" +Xh7+ ")+0.0007035052/" +Xh7+ ")\n"));
            }else if(0.25 < Xh7 && Xh7 <= 0.5){
                document.add(new Paragraph("V     = 1.56323-1.80696*log(Xg)-1.33458/Xh+0.276415*pow(log(Xg),2)+0.417135/pow(Xh,2)+1.39511*log(Xg)/Xh+0.0137129*pow(log(Xg),3)+0.0943597/pow(log(Xh),3)-0.402096*log(Xg)/pow(Xh,2)-0.101619*pow(log(Xg),2)/Xh\n"));
                document.add(new Paragraph("= 1.56323-1.80696*log(" +Xg7+ ")-1.33458/" +Xh7+ "+0.276415*pow(log(" +Xg7+ "),2)+0.417135/pow(" +Xh7+ ",2)+1.39511*log(" +Xg7+ ")/" +Xh7+ "+0.0137129*pow(log(" +Xg7+ "),3)+0.0943597/pow(log(" +Xh7+ "),3)-0.402096*log(" +Xg7+ ")/pow(" +Xh7+ ",2)-0.101619*pow(log(" +Xg7+ "),2)/" +Xh7+ "\n"));
            }else if(0.5 < Xh7 && Xh7 <= 1){
                document.add(new Paragraph("V     = -0.0213643-0.0763597/Xg+0.10299/Xh+0.725776/pow(Xg,2)-0.160603/pow(Xh,2)-0.0918061/(Xg*Xh)+0.472277/pow(Xg,3)+0.087353/pow(Xh,3)+0.527487/(Xg*pow(Xh,2))-0.980209/(pow(Xg,2)*Xh)\n"));
                document.add(new Paragraph("= -0.0213643-0.0763597/" +Xg7+ "+0.10299/" +Xh7+ "+0.725776/pow(" +Xg7+ ",2)-0.160603/pow(" +Xh7+ ",2)-0.0918061/(" +Xg7+ "*" +Xh7+ ")+0.472277/pow(" +Xg7+ ",3)+0.087353/pow(" +Xh7+ ",3)+0.527487/(" +Xg7+ "*pow(" +Xh7+ ",2))-0.980209/(pow(" +Xg7+ ",2)*" +Xh7+ ")\n"));
            }else{
                document.add(new Paragraph("V     = 0.00796687-0.220518/Xg+0.0602652/Xh+0.619818/pow(Xg,2)-0.223212/pow(Xh,2)+0.42192/(Xg*Xh)+0.0950195/pow(Xg,3)+0.209813/pow(Xh,3)-0.158821/(Xg*pow(Xh,2))-0.242056/(pow(Xg,2)*Xh)\n"));
                document.add(new Paragraph("= 0.00796687-0.220518/" +Xg7+ "+0.0602652/" +Xh7+ "+0.619818/pow(" +Xg7+ ",2)-0.223212/pow(" +Xh7+ ",2)+0.42192/(" +Xg7+ "*" +Xh7+ ")+0.0950195/pow(" +Xg7+ ",3)+0.209813/pow(" +Xh7+ ",3)-0.158821/(" +Xg7+ "*pow(" +Xh7+ ",2))-0.242056/(pow(" +Xg7+ ",2)*" +Xh7+ ")\n")); } }
        document.add(new Paragraph("= " +V7+ "  \n\n"));


        document.add(new Paragraph("e     = F/ho\n"));
        document.add(new Paragraph("= " +F7+ "/" +ho7+ "\n"));
        document.add(new Paragraph("= " +e7+ "  \n\n"));


        document.add(new Paragraph("d     = U*pow(go,2)*ho/V\n"));
        document.add(new Paragraph("= " +U7+ "*pow(" +go7+ ",2)*" +ho7+ "/" +V7+ "\n"));
        document.add(new Paragraph("= " +d7+ "  \n\n"));


        document.add(new Paragraph("L     = (t*e+1)/T+pow(t,3)/d\n"));
        document.add(new Paragraph("= (" +t7+ "*" +e7+ "+1)/" +T7+ "+pow(" +t7+ ",3)/" +d7+ "\n"));
        document.add(new Paragraph("= " +L7+ "  \n\n"));


        if(FT7 == 0){
            document.add(new Paragraph("f     = max(1,(0.0927779-0.0336633*Xg+0.964176*pow(Xg,2)+0.0566286*Xh+0.347074*pow(Xh,2)-4.18699*pow(Xh,3))/(1-0.00596093*Xg+1.62904*Xh+3.49329*pow(Xh,2)+1.39052*pow(Xh,3)))\n"));
            document.add(new Paragraph("= max(1,(0.0927779-0.0336633*" +Xg7+ "+0.964176*pow(" +Xg7+ ",2)+0.0566286*" +Xh7+ "+0.347074*pow(" +Xh7+ ",2)-4.18699*pow(" +Xh7+ ",3))/(1-0.00596093*" +Xg7+ "+1.62904*" +Xh7+ "+3.49329*pow(" +Xh7+ ",2)+1.39052*pow(" +Xh7+ ",3)))\n"));
            document.add(new Paragraph("= max(1," +(0.0927779-0.0336633*Xg7+0.964176*Math.pow(Xg7,2)+0.0566286*Xh7+0.347074*Math.pow(Xh7,2)-4.18699*Math.pow(Xh7,3))/(1-0.00596093*Xg7+1.62904*Xh7+3.49329*Math.pow(Xh7,2)+1.39052*Math.pow(Xh7,3))+ ")\n"));
            document.add(new Paragraph("= " +f7+ "  \n\n"));
        }else{
            document.add(new Paragraph("f     = " +f7+ "  \n\n"));
        }



        document.add(new Paragraph("H     = 0.785*pow(G,2)*P\n"));
        document.add(new Paragraph("= 0.785*pow(" +G7+ ",2)*" +P6+ "\n"));
        document.add(new Paragraph("= " +H7+ "  \n\n"));


        document.add(new Paragraph("HD    = 0.785*pow(B,2)*P\n"));
        document.add(new Paragraph("= 0.785*pow(" +B7+ ",2)*" +P6+ "\n"));
        document.add(new Paragraph("= " +HD7+ "  \n\n"));


        document.add(new Paragraph("HT    = H-HD\n"));
        document.add(new Paragraph("= " +H7+ "-" +HD7+ "\n"));
        document.add(new Paragraph("= " +HT7+ "  \n\n"));


        document.add(new Paragraph("HG    = Wo-H\n"));
        document.add(new Paragraph("= " +Wo7+ "-" +H7+ "\n"));
        document.add(new Paragraph("= " +HG7+ "  \n\n"));


        if(FT7 == 0 && PT7 == 0){
            document.add(new Paragraph("hD    = (C-B-g1)/2\n"));
            document.add(new Paragraph("= (" +C7+ "-" +B7+ "-" +g17+ ")/2\n"));
        }else if(FT7 == 0 && PT7 == 1){
            document.add(new Paragraph("hD    = (C-B+g1-2*go)/2\n"));
            document.add(new Paragraph("= (" +C7+ "-" +B7+ "+" +g17+ "-2*" +go7+ ")/2\n"));
        }else{
            document.add(new Paragraph("hD    = (C-B)/2\n"));
            document.add(new Paragraph("= (" +C7+ "-" +B7+ ")/2\n")); }
        document.add(new Paragraph("= " +hD7+ "  \n\n"));


        document.add(new Paragraph("hG    = (C-G)/2\n"));
        document.add(new Paragraph("= (" +C7+ "-" +G7+ ")/2\n"));
        document.add(new Paragraph("= " +hG7+ "  \n\n"));


        document.add(new Paragraph("hT    = 0.5*(C-(B+G)/2\n"));
        document.add(new Paragraph("= 0.5*(" +C7+ "-(" +B7+ "+" +G7+ ")/2)\n"));
        document.add(new Paragraph("= " +hT7+ "  \n\n"));


        document.add(new Paragraph("Bs    = C/nb\n"));
        document.add(new Paragraph("= " +C7+ "/" +nb6+ "\n"));
        document.add(new Paragraph("= " +Bs7+ "  \n\n"));


        document.add(new Paragraph("Bsc   = max(1,sqrt(Bs/(2*db+t)))\n"));
        document.add(new Paragraph("= max(1,sqrt(" +Bs7+ "/(2*" +db7+ "+" +t7+ ")))\n"));
        document.add(new Paragraph("= max(1," +(Math.sqrt(Bs7/(2*db7+t7)))+ ")\n"));
        document.add(new Paragraph("= " +Bsc7+ "  \n\n"));


        document.add(new Paragraph("I     = 0.0874*L*pow(go,2)*ho*B/V\n"));
        document.add(new Paragraph("= 0.0874*" +L7+ "*pow(" +go7+ ",2)*" +ho7+ "*" +B7+ "/" +V7+ "\n"));
        document.add(new Paragraph("= " +I7+ "  \n\n"));


        document.add(new Paragraph("AR    = 0.5*(A-B)\n"));
        document.add(new Paragraph("= 0.5*(" +A7+ "-" +B7+ ")\n"));
        document.add(new Paragraph("= " +AR7+ "  \n\n"));


        document.add(new Paragraph("Gavg  = 0.5*(go+g1)\n"));
        document.add(new Paragraph("= 0.5*(" +go7+ "+" +g17+ ")\n"));
        document.add(new Paragraph("= " +Gavg7+ "  \n\n"));


        if(t7 >= Gavg7){
            document.add(new Paragraph("AA    = AR = " +AA7+ "\n"));
            document.add(new Paragraph("BB    = t = " +BB7+ "\n"));
            document.add(new Paragraph("CC    = h = " +CC7+ "\n"));
            document.add(new Paragraph("DDG   = Gavg = " +DDG7+ "  \n\n"));
        }else{
            document.add(new Paragraph("AA    = h+t = " +h7+ "+" +t7+ " = " +AA7+ "\n"));
            document.add(new Paragraph("BB    = Gavg = " +BB7+ "\n"));
            document.add(new Paragraph("CC    = AR-Gavg = " +AR7+ "-" +Gavg7+ " = " +CC7+ "\n"));
            document.add(new Paragraph("DDG   = t = " +DDG7+ "  \n\n")); }


        document.add(new Paragraph("KAB   = (AA*pow(BB,3))*(0.33333-0.21*BB/AA*(1-0.08333*pow(BB/AA,4)))\n"));
        document.add(new Paragraph("= (" +AA7+ "*pow(" +BB7+ ",3))*(0.33333-0.21*" +BB7+ "/" +AA7+ "*(1-0.08333*pow(" +BB7+ "/" +AA7+ ",4)))\n"));
        document.add(new Paragraph("= " +KAB7+ "  \n\n"));


        document.add(new Paragraph("KCD   = (CC*pow(DDG,3))*(0.33333-(0.105*(DDG/CC)*(1-0.005208*pow((DDG/CC),4))))\n"));
        document.add(new Paragraph("= (" +CC7+ "*pow(" +DDG7+ ",3))*(0.33333-(0.105*(" +DDG7+ "/" +CC7+ ")*(1-0.005208*pow((" +DDG7+ "/" +CC7+ "),4))))\n"));
        document.add(new Paragraph("= " +KCD7+ "  \n\n"));


        document.add(new Paragraph("Ip    = KAB+KCD\n"));
        document.add(new Paragraph("= " +KAB7+ "+" +KCD7+ "\n"));
        document.add(new Paragraph("= " +Ip7+ "  \n\n"));


        document.add(new Paragraph("Moe   = 4*ME*(I/(0.3846*Ip+1))*(hD/(C-2*hD))+FA*hD\n"));
        document.add(new Paragraph("= 4*" +ME7+ "*(" +I7+ "/(0.3846*" +Ip7+ "+1))*(" +hD7+ "/(" +C7+ "-2*" +hD7+ "))+" +FA7+ "*" +hD7+ "\n"));
        document.add(new Paragraph("= " +Moe7+ "  \n\n"));


        if( PT7 == 0 ){
            document.add(new Paragraph("Mo    = abs(((HD*hD+HT*hT+HG*hG)*Bsc+Moe)*Fs)\n"));
            document.add(new Paragraph("= abs(((" +HD7+ "*" +hD7+ "+" +HT7+ "*" +hT7+ "+" +HG7+ "*" +hG7+ ")*" +Bsc7+ "+" +Moe7+ ")*1)\n"));
            document.add(new Paragraph("= " +Mo7+ "  \n\n"));

            document.add(new Paragraph("Mg    = Wg*(C-G)*Bsc*Fs/2\n"));
            document.add(new Paragraph("= " +Wg7+ "*(" +C7+ "-" +G7+ ")*" +Bsc7+ "*1/2\n"));
        }else{
            document.add(new Paragraph("Mo    = abs((HD*(hD-hG)+HT*(hT-hG)+Moe)*Fs)\n"));
            document.add(new Paragraph("= abs((" +HD7+ "*(" +hD7+ "-" +hG7+ ")+" +HT7+ "*(" +hT7+ "-" +hG7+ ")+" +Moe7+ ")*1)\n"));
            document.add(new Paragraph("= " +Mo7+ "  \n\n"));

            document.add(new Paragraph("Mg    = Wg*hG*Fs\n"));
            document.add(new Paragraph("= " +Wg7+ "*" +hG7+ "*1\n")); }
        document.add(new Paragraph("= " +Mg7+ "  \n\n"));


        document.add(new Paragraph("SHo   = f*Mo/(L*pow(g1,2)*B)\n"));
        document.add(new Paragraph("= " +f7+ "*" +Mo7+ "/(" +L7+ "*pow(" +g17+ ",2)*" +B7+ ")\n"));
        document.add(new Paragraph("= " +SHo7+ "  \n\n"));


        document.add(new Paragraph("SHg   = f*Mo/(L*pow(g1,2)*B)\n"));
        document.add(new Paragraph("= " +f7+ "*" +Mg7+ "/(" +L7+ "*pow(" +g17+ ",2)*" +B7+ ")\n"));
        document.add(new Paragraph("= " +SHg7+ "  \n\n"));


        document.add(new Paragraph("SRo   = (1.33*t*e+1)*Mo/(L*pow(t,2)*B)\n"));
        document.add(new Paragraph("= (1.33*" +t7+ "*" +e7+ "+1)*" +Mo7+ "/(" +L7+ "*pow(" +t7+ ",2)*" +B7+ ")\n"));
        document.add(new Paragraph("= " +SRo7+ "  \n\n"));


        document.add(new Paragraph("SRg   = (1.33*t*e+1)*Mo/(L*pow(t,2)*B)\n"));
        document.add(new Paragraph("= (1.33*" +t7+ "*" +e7+ "+1)*" +Mg7+ "/(" +L7+ "*pow(" +t7+ ",2)*" +B7+ ")\n"));
        document.add(new Paragraph("= " +SRg7+ "  \n\n"));


        if( PT7 == 0 ){
            document.add(new Paragraph("STo    = Y*Mo/(pow(t,2)*B)-Z*SRo\n"));
            document.add(new Paragraph("= " +Y7+ "*" +Mo7+ "/(pow(" +t7+ ",2)*" +B7+ ")-" +Z7+ "*" +SRo7+ "\n"));
            document.add(new Paragraph("= " +STo17+ "  \n\n"));

            document.add(new Paragraph("STg    = Y*Mg/(pow(t,2)*B)-Z*SRg\n"));
            document.add(new Paragraph("= " +Y7+ "*" +Mg7+ "/(pow(" +t7+ ",2)*" +B7+ ")-" +Z7+ "*" +SRg7+ "\n"));
            document.add(new Paragraph("= " +STg17+ "  \n\n"));
        }else{
            document.add(new Paragraph("STo1    = Y*Mo/(pow(t,2)*B)-Z*SRo*(0.67*t*e+1)/(1.33*t*e+1)\n"));
            document.add(new Paragraph("= " +Y7+ "*" +Mo7+ "/(pow(" +t7+ ",2)*" +B7+ ")-" +Z7+ "*" +SRo7+ "*(0.67*" +t7+ "*" +e7+ "+1)/(1.33*" +t7+ "*" +e7+ "+1)\n"));
            document.add(new Paragraph("= " +STo17+ "  \n\n"));

            document.add(new Paragraph("STg1    = Y*Mg/(pow(t,2)*B)-Z*SRg*(0.67*t*e+1)/(1.33*t*e+1)\n"));
            document.add(new Paragraph("= " +Y7+ "*" +Mg7+ "/(pow(" +t7+ ",2)*" +B7+ ")-" +Z7+ "*" +SRg7+ "*(0.67*" +t7+ "*" +e7+ "+1)/(1.33*" +t7+ "*" +e7+ "+1)\n"));
            document.add(new Paragraph("= " +STg17+ "  \n\n"));

            document.add(new Paragraph("STo2    = Mo/(pow(t,2)*B)*(Y-2*pow(K,2)*(0.67*t*e+1)/((pow(K,2)-1)*L))\n"));
            document.add(new Paragraph("= " +Mo7+ "/(pow(" +t7+ ",2)*" +B7+ ")*(" +Y7+ "-2*pow(" +K7+ ",2)*(0.67*" +t7+ "*" +e7+ "+1)/((pow(" +K7+ ",2)-1)*" +L7+ "))\n"));
            document.add(new Paragraph("= " +STo27+ "  \n\n"));

            document.add(new Paragraph("STg2    = Mg/(pow(t,2)*B)*(Y-2*pow(K,2)*(0.67*t*e+1)/((pow(K,2)-1)*L))\n"));
            document.add(new Paragraph("= " +Mg7+ "/(pow(" +t7+ ",2)*" +B7+ ")*(" +Y7+ "-2*pow(" +K7+ ",2)*(0.67*" +t7+ "*" +e7+ "+1)/((pow(" +K7+ ",2)-1)*" +L7+ "))\n"));
            document.add(new Paragraph("= " +STg27+ "  \n\n")); }




        //Flange stress acceptance criteria
        //For operating condition
        document.newPage();
        document.add(new Paragraph("Flange Calculator - 2019  ©Sushil Champ"));
        document.add(new Paragraph("Report created on: "+new Date()));
        document.add(new Paragraph("****************************************************************************************************************\n"));


        document.add(new Paragraph("Checking flange stress acceptance criterion for operating condition:\n\n"));
        if(PT7 == 0){               ///////For integral and loose type flange
            //First condition
            if(FT7 == 0){
                document.add(new Paragraph("       SHo <= min(1.5*Sfo,2.5*Sno)\n"));
                document.add(new Paragraph("==>  " +SHo7+ " <= min(1.5*" +Sfo6+ ",2.5*" +Sno6+ ")\n"));
                document.add(new Paragraph("==>  " +SHo7+ " <= min(" +(1.5*Sfo6)+ "," +(2.5*Sno6)+ ")\n"));
                document.add(new Paragraph("==>  " +SHo7+ " <= " +Math.min(1.5*Sfo6,2.5*Sno6)+ "\n"));
                if(SHo7 <= Math.min(1.5*Sfo6,2.5*Sno6)){
                    document.add(new Paragraph("==>  Satisfied\n"));
                }else{
                    document.add(new Paragraph("==>  Not Satisfied\n")); }
            }else{
                document.add(new Paragraph("      SHo <= 1.5*Sfo\n"));
                document.add(new Paragraph("==>  " +SHo7+ " <= 1.5*" +Sfo6+ "\n"));
                document.add(new Paragraph("==>  " +SHo7+ " <= " +(1.5*Sfo6)+ "\n"));
                if(SHo7 <= 1.5*Sfo6){
                    document.add(new Paragraph("==>  Satisfied\n"));
                }else{
                    document.add(new Paragraph("==>  Not Satisfied\n")); } }
            document.add(new Paragraph("\n"));

            //Second condition
            document.add(new Paragraph("       SRo <= Sfo\n"));
            document.add(new Paragraph("==>  " +SRo7+ " <= " +Sfo6+ "\n"));
            if(SRo7 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Third condition
            document.add(new Paragraph("       STo <= Sfo\n"));
            document.add(new Paragraph("==>  " +STo17+ " <= " +Sfo6+ "\n"));
            if(STo17 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Fourth condition
            document.add(new Paragraph("       ((SHo+SRo)/2) <= Sfo\n"));
            document.add(new Paragraph("==>  ((" +SHo7+ "+" +SRo7+ ")/2) <= " +Sfo6+ "\n"));
            document.add(new Paragraph("==>  " +(SHo7+SRo7)/2+ " <= " +Sfo6+ "\n"));
            if((SHo7+SRo7)/2 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Fifth condition
            document.add(new Paragraph("       ((SHo+STo)/2) <= Sfo\n"));
            document.add(new Paragraph("==>  ((" +SHo7+ "+" +STo17+ ")/2) <= " +Sfo6+ "\n"));
            document.add(new Paragraph("==>  " +(SHo7+STo17)/2+ " <= " +Sfo6+ "\n"));
            if((SHo7+STo17)/2 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

        }else{                  /////////For reverse integral and loose type flange
            //First condition
            document.add(new Paragraph("       SHo <= 1.5*Sfo\n"));
            document.add(new Paragraph("==>  " +SHo7+ " <= 1.5*" +Sfo6+ "\n"));
            document.add(new Paragraph("==>  " +SHo7+ " <= " +(1.5*Sfo6)+ "\n"));
            if(SHo7 <= 1.5*Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Second condition
            document.add(new Paragraph("       SRo <= Sfo\n"));
            document.add(new Paragraph("==>  " +SRo7+ " <= " +Sfo6+ "\n"));
            if(SRo7 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Third condition
            document.add(new Paragraph("       STo <= Sfo\n"));
            document.add(new Paragraph("==>  " +STo17+ " <= " +Sfo6+ "\n"));
            if(STo17 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Fourth condition
            document.add(new Paragraph("       ((SHo+SRo)/2) <= Sfo\n"));
            document.add(new Paragraph("==>  ((" +SHo7+ "+" +SRo7+ ")/2) <= " +Sfo6+ "\n"));
            document.add(new Paragraph("==>  " +(SHo7+SRo7)/2+ " <= " +Sfo6+ "\n"));
            if((SHo7+SRo7)/2 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Fifth condition
            document.add(new Paragraph("       ((SHo+STo)/2) <= Sfo\n"));
            document.add(new Paragraph("==>  ((" +SHo7+ "+" +STo17+ ")/2) <= " +Sfo6+ "\n"));
            document.add(new Paragraph("==>  " +(SHo7+STo17)/2+ " <= " +Sfo6+ "\n"));
            if((SHo7+STo17)/2 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Sixth condition
            document.add(new Paragraph("       STo <= Sfo\n"));
            document.add(new Paragraph("==>  " +STo27+ " <= " +Sfo6+ "\n"));
            if(STo27 <= Sfo6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));
        }



        //Flange stress acceptance criteria
        //For gasket seating condition
        document.newPage();
        document.add(new Paragraph("Flange Calculator - 2019  ©Sushil Champ"));
        document.add(new Paragraph("Report created on: "+new Date()));
        document.add(new Paragraph("****************************************************************************************************************\n"));


        document.add(new Paragraph("Checking flange stress acceptance criterion for gasket seating condition:\n\n"));
        if(PT7 == 0){               ///////For integral and loose type flange
            //First condition
            if(FT7 == 0){
                document.add(new Paragraph("       SHg <= min(1.5*Sfg,2.5*Sng)\n"));
                document.add(new Paragraph("==>  " +SHg7+ " <= min(1.5*" +Sfg6+ ",2.5*" +Sng6+ ")\n"));
                document.add(new Paragraph("==>  " +SHg7+ " <= min(" +(1.5*Sfg6)+ "," +(2.5*Sng6)+ ")\n"));
                document.add(new Paragraph("==>  " +SHg7+ " <= " +Math.min(1.5*Sfg6,2.5*Sng6)+ "\n"));
                if(SHg7 <= Math.min(1.5*Sfg6,2.5*Sng6)){
                    document.add(new Paragraph("==>  Satisfied\n"));
                }else{
                    document.add(new Paragraph("==>  Not Satisfied\n")); }
            }else{
                document.add(new Paragraph("      SHg <= 1.5*Sfg\n"));
                document.add(new Paragraph("==>  " +SHg7+ " <= 1.5*" +Sfg6+ "\n"));
                document.add(new Paragraph("==>  " +SHg7+ " <= " +(1.5*Sfg6)+ "\n"));
                if(SHg7 <= 1.5*Sfg6){
                    document.add(new Paragraph("==>  Satisfied\n"));
                }else{
                    document.add(new Paragraph("==>  Not Satisfied\n")); } }
            document.add(new Paragraph("\n"));

            //Second condition
            document.add(new Paragraph("       SRg <= Sfg\n"));
            document.add(new Paragraph("==>  " +SRg7+ " <= " +Sfg6+ "\n"));
            if(SRg7 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Third condition
            document.add(new Paragraph("       STg <= Sfg\n"));
            document.add(new Paragraph("==>  " +STg17+ " <= " +Sfg6+ "\n"));
            if(STg17 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Fourth condition
            document.add(new Paragraph("       ((SHg+SRg)/2) <= Sfg\n"));
            document.add(new Paragraph("==>  ((" +SHg7+ "+" +SRg7+ ")/2) <= " +Sfg6+ "\n"));
            document.add(new Paragraph("==>  " +(SHg7+SRg7)/2+ " <= " +Sfg6+ "\n"));
            if((SHg7+SRg7)/2 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Fifth condition
            document.add(new Paragraph("       ((SHg+STg)/2) <= Sfg\n"));
            document.add(new Paragraph("==>  ((" +SHg7+ "+" +STg17+ ")/2) <= " +Sfg6+ "\n"));
            document.add(new Paragraph("==>  " +(SHg7+STg17)/2+ " <= " +Sfg6+ "\n"));
            if((SHg7+STg17)/2 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

        }else{                  /////////For reverse integral and loose type flange
            //First condition
            document.add(new Paragraph("       SHg <= 1.5*Sfg\n"));
            document.add(new Paragraph("==>  " +SHg7+ " <= 1.5*" +Sfg6+ "\n"));
            document.add(new Paragraph("==>  " +SHg7+ " <= " +(1.5*Sfg6)+ "\n"));
            if(SHg7 <= 1.5*Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Second condition
            document.add(new Paragraph("       SRg <= Sfg\n"));
            document.add(new Paragraph("==>  " +SRg7+ " <= " +Sfg6+ "\n"));
            if(SRg7 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Third condition
            document.add(new Paragraph("       STg <= Sfg\n"));
            document.add(new Paragraph("==>  " +STg17+ " <= " +Sfg6+ "\n"));
            if(STg17 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Fourth condition
            document.add(new Paragraph("       ((SHg+SRg)/2) <= Sfg\n"));
            document.add(new Paragraph("==>  ((" +SHg7+ "+" +SRg7+ ")/2) <= " +Sfg6+ "\n"));
            document.add(new Paragraph("==>  " +(SHg7+SRg7)/2+ " <= " +Sfg6+ "\n"));
            if((SHg7+SRg7)/2 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Fifth condition
            document.add(new Paragraph("       ((SHg+STg)/2) <= Sfg\n"));
            document.add(new Paragraph("==>  ((" +SHg7+ "+" +STg17+ ")/2) <= " +Sfg6+ "\n"));
            document.add(new Paragraph("==>  " +(SHg7+STg17)/2+ " <= " +Sfg6+ "\n"));
            if((SHg7+STg17)/2 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));

            //Sixth condition
            document.add(new Paragraph("       STg <= Sfg\n"));
            document.add(new Paragraph("==>  " +STg27+ " <= " +Sfg6+ "\n"));
            if(STg27 <= Sfg6){
                document.add(new Paragraph("==>  Satisfied\n"));
            }else{
                document.add(new Paragraph("==>  Not Satisfied\n")); }
            document.add(new Paragraph("\n"));
        }




        //Flange rigidity criteria
        //For both operating and gasket seating condition
        document.newPage();
        document.add(new Paragraph("Flange Calculator - 2019  ©Sushil Champ"));
        document.add(new Paragraph("Report created on: "+new Date()));
        document.add(new Paragraph("****************************************************************************************************************\n"));


        document.add(new Paragraph("Checking flange rigidity criteria for operating condition:\n"));
        document.add(new Paragraph("Jo   = 52.14*V*Mo/(L*Eyo*pow(go,2)*KR*ho) <= 1\n"));
        document.add(new Paragraph("==>  52.14*" +V7+ "*" +Mo7+ "/(" +L7+ "*" +Eyo7+ "*pow(" +go7+ ",2)*" +KR7+ "*" +ho7+ ") <= 1\n"));
        document.add(new Paragraph("==>  " +Jo7+ " <= 1\n"));
        if(Jo7 <= 1){
            document.add(new Paragraph("==>  Satisfied\n"));
        }else{
            document.add(new Paragraph("==>  Not Satisfied\n")); }
        document.add(new Paragraph("\n"));

        document.add(new Paragraph("Checking flange rigidity criteria for gasket seating condition:\n"));
        document.add(new Paragraph("Jg   = 52.14*V*Mg/(L*Eyg*pow(go,2)*KR*ho) <= 1\n"));
        document.add(new Paragraph("==>  52.14*" +V7+ "*" +Mg7+ "/(" +L7+ "*" +Eyg7+ "*pow(" +go7+ ",2)*" +KR7+ "*" +ho7+ ") <= 1\n"));
        document.add(new Paragraph("==>  " +Jg7+ " <= 1\n"));
        if(Jg7 <= 1){
            document.add(new Paragraph("==>  Satisfied\n"));
        }else{
            document.add(new Paragraph("==>  Not Satisfied\n")); }
        document.add(new Paragraph("\n"));






        document.newPage();
        document.add(new Paragraph("Flange Calculator - 2019  ©Sushil Champ"));
        document.add(new Paragraph("Report created on: "+new Date()));
        document.add(new Paragraph("****************************************************************************************************************\n"));
        document.add(new Paragraph("3.   Results in Summary:\n\n"));
        document.add(new Paragraph("All the results of flange calculations are summarized as follows:\n\n"));
        document.add(ptresult);
        document.add(new Paragraph("\nOC: Operating Condition\n"));
        document.add(new Paragraph("GSC: Gasket Seating Condition\n"));
        document.add(new Paragraph("FCS: Flange Cross-Section\n\n"));
        if(areacheck == 1 && Jocheck == 1 && Jgcheck == 1 && stressOC == 1 && stressGSC == 1){
            document.add(new Paragraph("All parameters are passed.\nThe Flange is safe for provided design conditions.\n\n"));
        }else{
            document.add(new Paragraph("One or more parameters for selected flange are not satisfied with provided design conditions.\nThe Flange is not safe for provided design conditions.\n\n"));
            if(areacheck == 2){
                document.add(new Paragraph("Increase bolt size / bolt number or reduce bolt loads.\n"));
            }else if(stressOC == 2 || stressGSC == 2){
                document.add(new Paragraph("It is recommended to increase flange strength (select flange material of higher strength).\n"));
            }else{
                document.add(new Paragraph("Flange has lower rigidity than required.\n"));
                document.add(new Paragraph("Following solutions are recommended:\n"));
                document.add(new Paragraph("     Select flange of larger size/class.\n"));
                document.add(new Paragraph("     Select flange material of higher strength.\n"));
                document.add(new Paragraph("     Increase small hub thickness.\n"));
                document.add(new Paragraph("     Increase hub length.\n"));
            }
        }



        document.close();      //closes file after finishing up creating file
        previewPdf();          //calls for the function to open this pdf
    }



    //This is function #6-----This function views pdf if there is pdf viewer installed in the system
    private void previewPdf() {
        showMessageOKCancel("Do you want to preview the PDF file?",new DialogInterface.OnClickListener() {
            @Override
            public void onClick(DialogInterface dialog, int which) {
                PackageManager packManager = getPackageManager();
                Intent testIntent = new Intent(Intent.ACTION_VIEW);
                testIntent.setType("application/pdf");
                java.util.List<ResolveInfo> lists = packManager.queryIntentActivities(testIntent,PackageManager.MATCH_DEFAULT_ONLY);
                if (lists.size() > 0) {
                    Intent intent = new Intent();
                    intent.setAction(Intent.ACTION_VIEW);
                    Uri uri = Uri.fromFile(pdfFile);
                    intent.setDataAndType(uri, "application/pdf");
                    startActivity(intent);
                }else{
                    Toast.makeText(MainActivity.this, "Download a PDF Viewer to see the generated PDF", Toast.LENGTH_SHORT).show();
                }
            }
        });
    }



}