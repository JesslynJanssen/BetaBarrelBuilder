import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import g4p_controls.*; 
import controlP5.*; 
import processing.opengl.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class BetaBarrelGUI extends PApplet {

/*
 This program demonstrates the the GSketchPad control. This control enables
 drawings to be displayed on panels.
 
 for Processing V2 and V3
 (c) 2015 Peter Lager
 
 */



GPanel pnl;
GSketchPad spad;
PGraphics pg;


ControlP5 cp5;
int TurnProtein = 8;
int Strand_Number = 14;
float Radius=14.87f;
int Residues=15; 
float Strand_Angle=1.82f; 

Slider abc;



String fileName = "ideal (1).pdb";

PrintWriter output; 

int strandNumber=14; //strand number 
int numdata=0;
float allLines=0;
float curve=1.8233f; 
float radius=14.87f; 
float shearLength=39.321f; 
float barrelcenterX=-33.197448f; 
float barreltopZ= -93.887f; 
float barrelcenterY=-32.738695f; 
float CAstrand=15; 
static float x=0; 
static float y=0; 
static float z=0; 
float px=0; 
float py=0; 
float pz=0; 
float lpx=0; 
float lpy=0; 
float lpz=0; 
float scale=3;  //scaling just for viewing the protein
boolean osc=false;
float dev_ang = 0 * (PI/180);
float dev_rad = radius*0.05f; 
boolean valid=false;
float rotationValue=0;  //in radians
float CAangle=.0019f;
boolean print_pdb = false;
float tilt_ang = 0;





public void setup() {
 
 //   output = createWriter("pdb.txt");    //create a file to export PDB data to
  // Now create and clear graphic to be used (JAVA parameter
  // parameter is only needed for Processing 1.5.1)
 pg = createGraphics(500,500, OPENGL);
  // Create the G4P control to position and display the graphic
 spad = new GSketchPad(this, 0, 0, 500,500);
  // Set the graphic for this control. 
  // The graphic will be scaled to fit the control.
 spad.setGraphic(pg);
  // Create the panel and add the controls
 pnl = new GPanel(this, 10, 10, 500, 100, "Sketch Test");
 pnl.addControl(spad, 0, 0);
  // Expand the panel
 pnl.setCollapsed(false);
  
  
  cp5 = new ControlP5(this);


    // create another slider with tick marks, now without
  // default value, the initial value will be set according to
  // the value of variable sliderTicks2 then.
  cp5.addSlider("Strand_Number")
     .setPosition(550,240)
     .setSize(20,100)
     .setRange(8,24)
     .setNumberOfTickMarks(9);
     
     cp5.addSlider("Residues")
     .setPosition(640,240)
     .setSize(20,100)
     .setRange(10,35)
     .setNumberOfTickMarks(26);
     
     
      cp5.addSlider("Radius")
     .setPosition(720,240)
     .setSize(20,100)
     .setRange(8,24);
    
      cp5.addSlider("Strand_Angle")
     .setPosition(800,240)
     .setSize(20,100)
     .setRange(.7f,2);
    
       cp5.addSlider("Tilt")
     .setPosition(510,10)
     .setSize(20,100)
     .setRange(0,200);
    
    
     // add a vertical slider
  cp5.addSlider("Turn_Protein")
     .setPosition(160,515)
     .setSize(200,20)
     .setRange(0,200)
     .setValue(128)
     ; 
     
     cp5.addButton("Export_PDB")
     .setValue(0)
     .setPosition(160,570)
     .setSize(200,19)
     ;
}

float[][] k={
{4.088f,  1.575f,  1.406f},
{5.053f,  1.193f,  -0.413f},  
{4.115f,  0.546f,  2.524f},  
{3.786f,  0.837f, 3.673f},  
{4.517f,  -0.688f,  2.171f},  
{4.673f,  -1.754f,  3.142f}  
};

public void Turn_Protein(float rotate) {
  rotationValue = rotate;
  
}

public void Strand_Number(int number) {
  strandNumber= number;
  
}



public void Radius(float value) {
  radius= value;
  
}

public void Strand_Angle(float angle) {
  curve=angle;
  
}

public void Tilt(float rotate) {
 tilt_ang = rotate;
}

public void Residues(int num) {
  CAstrand= num;
  CAangle=curve/num;
 float xCAlength= sin(CAangle)*3.5f;
 float zLength= sqrt(((3.5f*3.5f)-(xCAlength*xCAlength)));            //pythagoras
 shearLength= num*zLength;
}

public void draw() {
 background(150,150,150);
   updateGraphic();
 if(print_pdb){
  output.flush(); //writes the remaining data to the file   
  output.close(); //Finishes the file 
  print_pdb=false;
 }
}



public void backbone_draw_p2p(float start_x, float start_y, float start_z, float end_x, float end_y, float end_z){

 if(start_x !=0){
  float[][] coef={
{0,  0,  0},
{0.046153846f,  0.309101832f,  0.644009217f},
{-0.516239316f,  0.221688195f,  1.305875576f},
{0.733333333f,  0.679783719f,  0.440668203f},
{1.165811966f,  0.735656353f,  -0.093894009f},
{1,  1,  1}};

    float Cx  =  coef[1][0] *(end_x-start_x) +start_x ;
    float Ox  =  coef[2][0] *(end_x-start_x) +start_x ;
    float Nx  =  coef[3][0] *(end_x-start_x) +start_x ;
    float Hx  =  coef[4][0] *(end_x-start_x) +start_x ;

    float Cy  =  coef[1][1] *(end_y-start_y) +start_y ;
    float Oy  =  coef[2][1] *(end_y-start_y) +start_y ;
    float Ny  =  coef[3][1] *(end_y-start_y) +start_y ;
    float Hy  =  coef[4][1] *(end_y-start_y) +start_y ;
    
    float Cz  =  coef[1][2] *(end_z-start_z) +start_z ;
    float Oz  =  coef[2][2] *(end_z-start_z) +start_z ;
    float Nz  =  coef[3][2] *(end_z-start_z) +start_z ;
    float Hz  =  coef[4][2] *(end_z-start_z) +start_z ;
    
    
    
    pg.strokeWeight(4); 
    pg.point(scale*start_x,scale*start_y,scale*start_z);
    pg.point(scale*Cx,scale*Cy,scale*Cz);
    pg.point(scale*Ox,scale*Oy,scale*Oz);
    pg.point(scale*Nx,scale*Ny,scale*Nz);
    pg.point(scale*Hx,scale*Hy,scale*Hz);
    pg.point(scale*end_x,scale*end_y,scale*end_z);
    
    pg.strokeWeight(1); 
    pg.line(scale*start_x,scale*start_y,scale*start_z,scale*Cx,scale*Cy,scale*Cz);
    pg.line(scale*Ox,scale*Oy,scale*Oz,scale*Cx,scale*Cy,scale*Cz);
    pg.line(scale*Cx,scale*Cy,scale*Cz,scale*Nx,scale*Ny,scale*Nz);
    pg.line(scale*Hx,scale*Hy,scale*Hz,scale*Nx,scale*Ny,scale*Nz);
    pg.line(scale*Nx,scale*Ny,scale*Nz,scale*end_x,scale*end_y,scale*end_z);
    
    if(print_pdb){
      String temp;
      
      temp ="ATOM"+"     "+"1"+" "+"CA"+"   "+"GLY"+" "+"A"+"  "+"2"+"     "+str(start_x) + "\t" + str(start_y) + "\t" + str(start_z);
      output.println(temp);
      
      temp ="ATOM"+"     "+"1"+" "+"C"+ "   "+"GLY"+" "+"A"+"  "+"2"+"     "+str(Cx) + "\t" + str(Cy) + "\t" + str(Cz);
      output.println(temp);
      
      temp = "ATOM"+"     "+"1"+" "+"O"+"   "+"GLY"+" "+"A"+"  "+"2"+"     "+str(Ox) + "\t" + str(Oy) + "\t" + str(Oz);
      output.println(temp);
      
      temp = "ATOM"+"     "+"1"+" "+"N"+"   "+"GLY"+" "+"A"+"  "+"2"+"     "+str(Nx) + "\t" + str(Ny) + "\t" + str(Nz);
      output.println(temp);
      
      temp = "ATOM"+"     "+"1"+" "+"H"+"   "+"GLY"+" "+"A"+"  "+"2"+"     "+str(Hx) + "\t" + str(Hy) + "\t" + str(Hz);
      output.println(temp);
      
     // pg.point(end_x,end_y,end_z);
      
    }
  }
}


public void loop_points(float start_x, float start_y, float start_z, float end_x, float end_y, float end_z, float z_offset){
  
  
  
  float point1_x=((end_x-start_x)/3)+start_x;
  float point1_y=((end_y-start_y)/3)+start_y+z_offset;
  float point1_z=((end_z-start_z)/3)+start_z;
  float point2_x=(2*(end_x-start_x)/3)+start_x;
  float point2_y=(2*(end_y-start_y)/3)+start_y+z_offset;
  float point2_z=(2*(end_z-start_z)/3)+start_z;
  
  backbone_draw_p2p(start_x,start_y,start_z, point1_x,point1_y,point1_z);
  backbone_draw_p2p(point1_x,point1_y,point1_z,point2_x,point2_y,point2_z);
  backbone_draw_p2p(point2_x,point2_y,point2_z,end_x,end_y,end_z);
  
}


  

// Add a line or ellipse to the sketchpad graphic
public void updateGraphic() {
  
  pg.beginDraw();
   
  // MAKE BETA BARREL MODEL
  pg.background(0);     
  pg.strokeWeight(4);  
  
  pg.translate(250,370,200);
 
 pg.rotateX(radians(-tilt_ang/3));
  pg.rotateY(radians(rotationValue));
  pg.translate(-scale*barrelcenterX,-scale*barrelcenterY,100);
  
  float psi=curve/CAstrand; 
  CAangle=psi;
  
  backbone_draw_p2p(0,0,0, 100,100,100);
  
  for(float phi=2; phi<=strandNumber; phi+=2){        
  pg.colorMode(HSB); //adding colour    
  pg.stroke((phi*256/strandNumber),255,255); //adding colour    
  
  x=barrelcenterX+((radius+dev_rad)*sin(dev_ang+PI+(curve+((2.0f*PI*(phi-1))/strandNumber))));   
  y=barrelcenterY+((radius+dev_rad)*cos(dev_ang+PI+(curve+((2.0f*PI*(phi-1))/strandNumber))));                  
  z=barreltopZ+((curve)*(shearLength/curve));
    
  loop_points(lpx,lpz,lpy,x,z,y, CAangle*(shearLength/curve));
  
  for(float theta=curve;theta>CAangle;theta-=CAangle){ 
    
   dev_ang = -1 * dev_ang;
   dev_rad = -1 * dev_rad;
  
    px=x;
    py=y;
    pz=z;
          
    x=barrelcenterX+((radius+dev_rad)*sin(dev_ang+PI+(theta+((2.0f*PI*(phi-1))/strandNumber))));   
    y=barrelcenterY+((radius+dev_rad)*cos(dev_ang+PI+(theta+((2.0f*PI*(phi-1))/strandNumber))));                  
    z=barreltopZ+((theta)*(shearLength/curve));
  //  point(scale*x,scale*z,scale*y);
    
    lpx=x;lpy=y;lpz=z;
    backbone_draw_p2p(px,pz,py, x,z,y);
    
    //output.println(z); //write the data to the file   

  }
  
   dev_ang = -1 * dev_ang;
   dev_rad = -1 * dev_rad;
   
   x=barrelcenterX+((radius+dev_rad)*sin(dev_ang+PI+(CAangle+((2.0f*PI*(phi))/strandNumber))));          
   y=barrelcenterY+((radius+dev_rad)*cos(dev_ang+PI+(CAangle+((2.0f*PI*(phi))/strandNumber))));        
   z=barreltopZ+((CAangle)*(shearLength/curve));       
  // point(scale*x,scale*z,scale*y); 
   
   loop_points(lpx,lpz,lpy,x,z,y, -CAangle*(shearLength/curve));
   
   
  for(float theta=CAangle;theta<curve;theta+=CAangle){   
    
    dev_ang = -1 * dev_ang;
    dev_rad = -1 * dev_rad;
    
    px=x;
    py=y;
    pz=z;
    
    x=barrelcenterX+((radius+dev_rad)*sin(dev_ang+PI+(theta+((2.0f*PI*(phi))/strandNumber))));          
    y=barrelcenterY+((radius+dev_rad)*cos(dev_ang+PI+(theta+((2.0f*PI*(phi))/strandNumber))));        
    z=barreltopZ+((theta)*(shearLength/curve));       
  //  point(scale*x,scale*z,scale*y);
    
    lpx=x;lpy=y;lpz=z;
    
    backbone_draw_p2p(px,pz,py, x,z,y);
    
  } 
   dev_ang = -1 * dev_ang;
   dev_rad = -1 * dev_rad;  
  } 
  pg.endDraw(); 
  
}

public void Export_PDB(int theValue) {
  println(theValue);
  output = createWriter("pdb.txt");    //create a file to export PDB data to
  print_pdb = true;
}
  public void settings() {  size(900, 700,OPENGL); }
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "BetaBarrelGUI" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
