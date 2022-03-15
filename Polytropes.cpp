//Generator of a number of random EOS
//check the speed of sound
//check for mass limit

using namespace std;
#include <cmath>
#include <vector>
#define  G 6.6732*TMath::Power(10,-8)
#define  c 2.9979*TMath::Power(10,10)
#define  Mo 1.987*TMath::Power(10,33)
#define  eo  TMath::Power(10,15)
#define  Mb  1.66*TMath::Power(10,-24)
#define  kappa c*c/(G*eo)





//needed for the crust
int crust(double p,vector<double> vp,vector<double> vd){//crust EOS
  int pos;
  for(int i=0;i<vp.size()-1;i++){
    if(vp.at(i)<=p && vp.at(i+1)>=p) return i;
  }
    
  }


///needed for the crust




double dmdr(double r,double m,double p,vector<double> vp,vector<double> ve){//energy density
  int pos;
  for(int i=0;i<vp.size()-1;i++){
    if(vp.at(i)<=p && vp.at(i+1)>=p) pos= i;
  }


  double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)));
  double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))*TMath::Log10(vp.at(pos+1)));
  double en=b*TMath::Power(p,a);
  return 4.*TMath::Pi()*r*r*en;
}


double dpdr(double r,double m,double p,vector<double> vp,vector<double> ve){
  int pos;
  for(int i=0;i<vp.size()-1;i++){
    if(vp.at(i)<=p && vp.at(i+1)>=p) pos= i;
  }

   double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)));
  double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))*TMath::Log10(vp.at(pos+1)));
  double en=b*TMath::Power(p,a);
  
  return -(en+p)*(m+4.*TMath::Pi()*r*r*r*p)/(r*(r-2.*m));

}//drho/dr

  
class findParameter {
public:
  double p0;
  double Gamma1;
  double Gamma2;
  double Gamma3;
  double K1;
  double K2;
  double K3;
  double rho1;
  double rho2;
  double rho0;
  double a1;
  double a2;
  double a3;

  
  void setParameters(double pa1,double pa2,double pa3,double pa4,double pa5, double par6, double par7){
    p0=TMath::Power(10,pa1)*kappa*G/(c*c*c*c);
    Gamma1=pa2;
    Gamma2=pa3;
    Gamma3=pa4;
    rho0=pa5*kappa*G/(c*c);
    rho1= par6*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    rho2= par7*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    K1=p0/TMath::Power(rho0,Gamma1);
    K2=K1*TMath::Power(rho1,Gamma1-Gamma2);
    K3=K2*TMath::Power(rho2,Gamma2-Gamma3);

  }

  void setAlphas(double ec, double dc,double pc){

    a1=-1-K1*TMath::Power(dc,Gamma1-1)/(Gamma1-1)+ec/dc;
    a2=-1-K2*TMath::Power(rho1,Gamma2-1)/(Gamma2-1)+(1+a1)+K1*TMath::Power(rho1,Gamma1-1)/(Gamma1-1);
    a3=-1-K3*TMath::Power(rho2,Gamma3-1)/(Gamma3-1)+(1+a2)+K2*TMath::Power(rho2,Gamma2-1)/(Gamma2-1);


  }


  //needed for the crust
  int crust(double p,vector<double> vp,vector<double> vd){//crust EOS
    int pos;
    for(int i=0;i<vp.size()-1;i++){
      if(vp.at(i)<=p && vp.at(i+1)>=p) return i;
    }
    
  }
  double dmdr(double r,double m,double p,vector<double> vp,vector<double> ve){//energy density
    int pos;
    for(int i=0;i<vp.size()-1;i++){
      if(vp.at(i)<=p && vp.at(i+1)>=p) pos= i;
    }


    double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)));
    double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))*TMath::Log10(vp.at(pos+1)));
    double en=b*TMath::Power(p,a);
    return 4.*TMath::Pi()*r*r*en;
  }


  double dpdr(double r,double m,double p,vector<double> vp,vector<double> ve){
    int pos;
    for(int i=0;i<vp.size()-1;i++){
      if(vp.at(i)<=p && vp.at(i+1)>=p) pos= i;
    }

    double a=(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)));
    double b=TMath::Power(10,TMath::Log10(ve.at(pos+1))-(TMath::Log10(ve.at(pos))-TMath::Log10(ve.at(pos+1)))/(TMath::Log10(vp.at(pos))-TMath::Log10(vp.at(pos+1)))*TMath::Log10(vp.at(pos+1)));
    double en=b*TMath::Power(p,a);
  
    return -(en+p)*(m+4.*TMath::Pi()*r*r*r*p)/(r*(r-2.*m));

  }
  //need for crust
  

  
  double pressure(double d){//EOS
    if(d<=rho1){
      return K1*TMath::Power(d,Gamma1);
    }
    else if(d>rho1 && d<rho2){
      return K2*TMath::Power(d,Gamma2);
    }
    else if(d>=rho2){
      return K3*TMath::Power(d,Gamma3);
    }
  }



  double h(double d){

    if(d<=rho1){
      return 1.+a1+K1*Gamma1*TMath::Power(d,Gamma1-1.)/(Gamma1-1.);
    }
    else if(d>rho1&&d<rho2){
      return 1.+a2+K2*Gamma2*TMath::Power(d,Gamma2-1.)/(Gamma2-1.);
    }
    else if(d>=rho2){
      return 1.+a3+K3*Gamma3*TMath::Power(d,Gamma3-1.)/(Gamma3-1.);
    }
  
 
  }
  double e(double d){//energy density

    if(d<=rho1){
      return (1+a1)*d+K1*TMath::Power(d,Gamma1)/(Gamma1-1.);
    }
    else if(d>rho1&&d<rho2){
      return (1+a2)*d+K2*TMath::Power(d,Gamma2)/(Gamma2-1.);
    }
    else if(d>=rho2){
      return (1+a3)*d+K3*TMath::Power(d,Gamma3)/(Gamma3-1.);
    }  
  }




  double f(double p,double en){
  
    return 1./(en*kappa*G/(c*c)+p*G*kappa/(c*c*c*c));
  }


  double enthalpy(vector<double> p,vector<double> en,double h0,int num,int data){

    double dx=(p.at(num)-p.at(data))/((num-data)*1.);
    double integr=0;
    double press=p.at(0);
    
    integr=f(p.at(data),en.at(data))+f(p.at(num),en.at(num));

    
    for(int i=1;i<(num-data)*1-1;i++){
      press+=dx;
      int pos=crust(press,p,en);

      double a=(TMath::Log10(en.at(pos))-TMath::Log10(en.at(pos+1)))/(TMath::Log10(p.at(pos))-TMath::Log10(p.at(pos+1)));
      double b=TMath::Power(10,TMath::Log10(en.at(pos+1))-(TMath::Log10(en.at(pos))-TMath::Log10(en.at(pos+1)))/(TMath::Log10(p.at(pos))-TMath::Log10(p.at(pos+1)))*TMath::Log10(p.at(pos+1)));
      double ener=b*TMath::Power(press,a);
      if(i%2==0){
	integr+=2.*f(press,ener);
      }
      else{
	integr+=4.*f(press,ener);
      }
    }

    integr=h0+integr*dx/(c*c*3.);
 
     
    return integr;
  }




  bool soundSpeed(double p,double d){

    if(d<=rho1){
      if(TMath::Sqrt(Gamma1*p/(e(d)+p))>1.){
	return false;
      }
      else{
	return true;
      }
    }
    else if(d>rho1&&d<rho2){
      if(TMath::Sqrt(Gamma2*p/(e(d)+p))>1.){
	return false;
      }
      else{
	return true;
      }
    }
    else if(d>=rho2){
      if(TMath::Sqrt(Gamma3*p/(e(d)+p))>1.){
	return false;
      }
      else{
	return true;
      }
    } 

  }


  
  double soundSpeedValue(double p,double d){

    if(d<=rho1){
      return TMath::Sqrt(Gamma1*p/(e(d)+p));
    }
    else if(d>rho1&&d<rho2){
      return TMath::Sqrt(Gamma2*p/(e(d)+p));
    }
    else if(d>=rho2){
      return TMath::Sqrt(Gamma3*p/(e(d)+p));
    } 

  }


  double maxMassOrRadius(vector<double> v_d,vector<double> v_p,vector<double> v_e,int m_or_r){

    int count=0;
    double stepsize=100./(TMath::Power(kappa,0.5));
    double r0=TMath::Power(10,-15)/(TMath::Power(kappa,0.5));
    double m0=TMath::Power(10,-15)*G/(c*c*TMath::Power(kappa,0.5));
    double d0=v_d.at(v_d.size()-1);
    double p0=v_p.at(v_p.size()-1);

    
    double r,m,d,p,k0,k1,k2,k3,l0,l1,l2,l3;

    d=d0;
    p=p0;
    m=m0;
  
    //RK4 Begining
    while(d>TMath::Power(10,0.9)*kappa*G/(c*c)){   
      if(count==0){
	/*
	  r=r0;//radius
	  m=m0;//mass
	  d=d0;//density
	  p=p0;//pressure
	*/
	r=r0;
	m=4*TMath::Pi()*d0*r0*r0*r0/3-TMath::Abs((TMath::Log10(d0)-TMath::Log10(v_d.size()-2))/(TMath::Log10(p0)-TMath::Log10(v_p.size()-2)))*8/15*3.14*3.14*d0/p0*(d0+p0)*(d0+3*p0)*TMath::Power(r0,5);
	d=d0-2/3*3.14*d0/p0*(d0+p0)*(d0+3*p0)*r0*r0;
	p=p0-TMath::Abs((TMath::Log10(d0)-TMath::Log10(v_d.at(v_p.size()-2)))/(TMath::Log10(p0)-TMath::Log10(v_p.at(v_p.size()-2))))*2*3.14*(d0+p0)*(d0+3*p0)*r0*r0/3;
      }
 
      else{
	
	double mold=m;
	double rold=r;
	double dold=d;
	double pold=p;
	double rnew,mnew,dnew,pnew;
	if(d<TMath::Power(10,14.5)*kappa*G/(c*c)) stepsize=10./(TMath::Power(kappa,0.5));
	if(d<TMath::Power(10,5.5)*kappa*G/(c*c)) stepsize=1./(TMath::Power(kappa,0.5));
	if(d<TMath::Power(10,3.5)*kappa*G/(c*c)) stepsize=0.01/(TMath::Power(kappa,0.5));
	if(d<TMath::Power(10,1.5)*kappa*G/(c*c)) stepsize=0.0001/(TMath::Power(kappa,0.5));
	
	k0=dmdr(rold,mold,pold,v_p,v_e);
	l0=dpdr(rold,mold,pold,v_p,v_e);

	if(pold+stepsize*0.5*l0<v_p.at(0))break;
	k1= dmdr(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,v_p,v_e);
	l1= dpdr(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,v_p,v_e);

	if(pold+stepsize*0.5*l1<v_p.at(0))break;
	k2 = dmdr(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,v_p,v_e);
	l2 = dpdr(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,v_p,v_e);

	if(pold+stepsize*l2<v_p.at(0))break;
	k3 = dmdr(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,v_p,v_e);
	l3 = dpdr(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,v_p,v_e);

	rnew=rold+stepsize;
	mnew=mold+(k0+2.*k1+2.*k2+k3)*stepsize/6.;
	pnew=pold+(l0+2.*l1+2.*l2+l3)*stepsize/6.;

	if(pnew<v_p.at(0))break;
	int pos=crust(pnew,v_p,v_d);

	double a=(TMath::Log10(v_d.at(pos))-TMath::Log10(v_d.at(pos+1)))/(TMath::Log10(v_p.at(pos))-TMath::Log10(v_p.at(pos+1)));
	double b=TMath::Power(10,TMath::Log10(v_d.at(pos+1))-(TMath::Log10(v_d.at(pos))-TMath::Log10(v_d.at(pos+1)))/(TMath::Log10(v_p.at(pos))-TMath::Log10(v_p.at(pos+1)))*TMath::Log10(v_p.at(pos+1)));
	dnew=b*TMath::Power(pnew,a);

 
      
	m=mnew;
	r=rnew;
	p=pnew;
	d=dnew;
	

      }
      count++;
    }

    if(m_or_r==0){
      return (m*c*c*TMath::Power(kappa,0.5))/(Mo*G);
    }else if(m_or_r==1){
      return (r*TMath::Power(kappa,0.5))/100000.;//Radius
    }
  }

  bool CheckDownRange(double p, double d){
    double r1= 2.5*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    double r2= 4.0*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    double Gam1=1.5;
    double Gam2=6;
    double Gam3=3;
    double Kappa1=p0/TMath::Power(rho0,Gam1);
    double Kappa2=Kappa1*TMath::Power(r1,Gam1-Gam2);
    double Kappa3=Kappa2*TMath::Power(r2,Gam2-Gam3);
    if(d<r1){
      if(p>Kappa1*TMath::Power(d,Gamma1)){
	return true;
      }
      else{
	return false;
      }
    }
    else if (d>r1 && d<r2){
      if(p>Kappa2*TMath::Power(d,Gamma2)){
	return true;
      }
      else{
	return false;
      }
    }
    else if (d>r2){
      if(p>Kappa3*TMath::Power(d,Gamma3)){
	return true;
      }
      else{
	return false;
      }
    }
	
  }
  
  bool CheckUpRange(double p, double d){
    double r1= 1.5*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    double r2= 2.0*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    double Gam1=4.5;
    double Gam2=5.5;
    double Gam3=3;
    double Kappa1=p0/TMath::Power(rho0,Gam1);
    double Kappa2=Kappa1*TMath::Power(r1,Gam1-Gam2);
    double Kappa3=Kappa2*TMath::Power(r2,Gam2-Gam3);
    if(d<r1){
       if(p<Kappa1*TMath::Power(d,Gamma1)){
	return true;
      }
      else{
	return false;
      }
    }
    else if (d>r1 && d<r2){
      if(p<Kappa2*TMath::Power(d,Gamma2)){
	return true;
      }
      else{
	return false;
      }
    }
    else if (d>r2){
      if(p<Kappa3*TMath::Power(d,Gamma3)){
	return true;
      }
      else{
	return false;
      }
    }
  }
  
};



void Polytropes(){//Main Program






  
  ofstream eosTable;
  eosTable.open("EosTable",ios::out);
  eosTable<<"EOS  log(p0) log(dyne/cm^2)  Gamma1  Gamma2  Gamma3 rho1 (g/cm^3) rho2 (g/cm^3)  Surface EOS"<<endl;
  bool Continue= false;

  //choose surface EOS
  vector<string> surfEOS;
  surfEOS.push_back("eosNV");
  surfEOS.push_back("eosFPS");
  surfEOS.push_back("eosBPS");
  int fileNum=0;



  //Number of random EOS
  int num;
  cout<<"Give the number of EOS"<<endl;
  cin>>num;
  
  //choose crust EOS
  int sEOSnum=0;//0 = NV   1 = FPS 2 = BPS

  //choose crust density
  double d_crust = 1.1*2.28*TMath::Power(10,14);//g/cm^3
  
  //maximum central mass density
  double Maxd0=8.3;//*2.28 10^14 g/cm^3

  //Maximum Mass
  double Max_mass_limit=1.97; //solar masses





  //Graphs
  vector<double> v1,v2,v3,v4;
  TCanvas* canvas =new TCanvas();
  canvas->cd();
  TImage *img = TImage::Create();

  //img->FromPad(5, 10, 10, 300, 200);
  img->FromPad(canvas);


  v3.push_back(6);
  v4.push_back(0);
  v3.push_back(6);
  v4.push_back(3);
  v3.push_back(18);
  v4.push_back(0);
  v3.push_back(18);
  v4.push_back(3);

  
  TGraph* m_r_cur =new TGraph(v4.size(),&v3.at(0),&v4.at(0));
  m_r_cur->SetTitle("Mass vs Radius Curves");
  m_r_cur->GetYaxis()->SetTitle("Mass (Mo)");
  m_r_cur->GetXaxis()->SetTitle("Radius (Km)");
  m_r_cur->SetLineColor(0);
  m_r_cur->SetMarkerColor(0);
  m_r_cur->SetMarkerStyle(20);
  m_r_cur->SetMarkerSize(1.5);
  
  m_r_cur->Draw("AP");

  

  
  //Random EOS Generation Begins
  while(fileNum<num){
       
    //Take the Surface EOS - Begining
    bool BoolSpeed=true;
    int count=0;
    vector<double> v_m,v_r,v_p,v_d,v_e,v_h;
    vector<double> vEnergyDensity,vEnthalpy,vDensity,vPressure,vSpeedofSound; //non dimensional
    ifstream myinputfile;
    myinputfile.open(surfEOS.at(sEOSnum),ios::in);
    string line;
    string a1,a2,a3,a4,a5;
    int data;
    while( getline(myinputfile,line) ){
	stringstream liness = stringstream(line);
	getline(liness,a1,' ');
	getline(liness,a2,' ');
	getline(liness,a3,' ');
	getline(liness,a4,' ');
	if(sEOSnum==0){getline(liness,a5,' ');}
	if((stod(a4)*1.66*TMath::Power(10,-24))>(d_crust)) break;
	v_e.push_back(stod(a1));
	v_p.push_back(stod(a2));
	v_h.push_back(stod(a3));
	v_d.push_back(stod(a4)*Mb);
	vPressure.push_back(stod(a2)*kappa*G/(c*c*c*c));
	vEnergyDensity.push_back(stod(a1)*kappa*G/(c*c));
	vEnthalpy.push_back(stod(a3)/(c*c));
	vDensity.push_back(stod(a4)*Mb*kappa*G/(c*c));
    }
    data=v_e.size();

    //Take the Surface EOS -END


    findParameter fp;
    double par1,par2,par3,par4,par5,par6,par7,d0;
    TRandom3* rand =new TRandom3();
    rand->SetSeed(0);
    par1=/*34.331;*/TMath::Log10(v_p.at(v_p.size()-1));//log(p0)
    par2=/*4.5;*/rand->Uniform(1,4.5);//Gamma1
    par3=/*2.835;*/rand->Uniform(0,8);//Gamma2
    par4=/*2.832;*/rand->Uniform(0.5,8);//Gamma3
    par5=v_d.at(v_d.size()-1)/*1.1*2.28*TMath::Power(10,14)*/;//rho0
    par6=/*2.198189621;*/rand->Uniform(1.5,Maxd0-1.0);//rho1*2.28 10^14 g/cm^3
    par7=/*4.385964912;*/rand->Uniform(par6,Maxd0-0.5);//rho2*2.28 10^14 g/cm^3
    d0=rand->Uniform(par7,Maxd0)/*Maxd0*/*2.28*TMath::Power(10,14)*kappa*G/(c*c);
    //call class
    fp.setParameters(par1,par2,par3,par4,par5,par6,par7);


    
    fp.setAlphas(vEnergyDensity.at(data-1),vDensity.at(data-1),vPressure.at(data-1));





    
    //EOS Begining
    double dens_stepsize=(d0-fp.rho0)/28.;
    double p0,d,p;
    d=1.1*2.28*TMath::Power(10,14)*(G*kappa)/(c*c)/*+fp.rho0*0.2*/;
    p=fp.pressure(d);
    bool CheckCont=true;
    bool up=true;
    bool down=true;
    
    while(d<=d0){
 
      v_d.push_back(d*c*c/(G*kappa));
      v_p.push_back(p*c*c*c*c/(G*kappa));
      v_e.push_back(fp.e(d)*c*c/(G*kappa));
	
      vPressure.push_back(p);
      vEnergyDensity.push_back(fp.e(d));
      vDensity.push_back(d);	
      //Finds N by using the politropic equation
      v_h.push_back(TMath::Log(fp.h(d))*(c*c));
      vEnthalpy.push_back(TMath::Log(fp.h(d)));
      vSpeedofSound.push_back(fp.soundSpeedValue(p,d));
      if( (v_h.at(v_h.size()-1)<v_h.at(v_h.size()-2))  || (v_p.at(v_p.size()-1)<v_p.at(v_p.size()-2)) || (v_d.at(v_d.size()-1)<v_d.at(v_d.size()-2)) || (v_e.at(v_e.size()-1)<v_e.at(v_e.size()-2))){
	CheckCont=false;
	break;
      }
	
      BoolSpeed=fp.soundSpeed(p,d);
      up=fp.CheckUpRange(p,d);
      down=fp.CheckDownRange(p,d);
      if(!BoolSpeed)break;
      if(!up)break;
      if(!down)break;
      d=d+dens_stepsize;
      p=fp.pressure(d);
      count++;
    }
    count=0;
    //EOS end


    

    //Check the speed of sound
    if(!up) continue;
    if(!down) continue;
    if(!BoolSpeed) continue;
    if(!CheckCont) continue;

    /*
    //Finds the N by integration
    
    for(int i=data;i<v_d.size();i++){
    double tempH=fp.enthalpy(v_p,v_e,v_h.at(data-1),i,data-1);
    v_h.push_back(tempH);
    vEnthalpy.push_back(tempH/(c*c));
    }
    */
    //Get e,p,h,d ENDS
   
    

    int stable_size=0;
    bool boolMax=true;
    double Maxmass=2.*Max_mass_limit;
    double mbefore,rbefore;


    //Find M-R curve and check if it passes the MaxMass-limit and Check for unstable part
    for(int dD=0;dD<vDensity.size()-data;dD++){
      double stepsize=100./(TMath::Power(kappa,0.5));
      double r0=TMath::Power(10,-15)/(TMath::Power(kappa,0.5));
      double m0=TMath::Power(10,-15)*G/(c*c*TMath::Power(kappa,0.5));
      double d0=vDensity.at(vDensity.size()-1-dD);//(TMath::Power(10,14.9)+dD*TMath::Power(10,15)*0.1)*kappa*G/(c*c);
      double p0=vPressure.at(vPressure.size()-1-dD);//pressure(d0);
      double maxSoundSpeed=0;

      
      //cout<<d0/(1.66*TMath::Power(10,-24)*kappa*G/(c*c))<<endl;

    
      double r,m,d,p,k0,k1,k2,k3,l0,l1,l2,l3;

      d=d0;
      p=p0;
      m=m0;
  
      //RK4 Begining
      while(d>TMath::Power(10,0.9)*kappa*G/(c*c)){   
	if(count==0){
	  r=r0;
	  m=4*TMath::Pi()*d0*r0*r0*r0/3-TMath::Abs((TMath::Log10(d0)-TMath::Log10(vDensity.size()-1-dD-1))/(TMath::Log10(p0)-TMath::Log10(vPressure.size()-1-dD-1)))*8/15*3.14*3.14*d0/p0*(d0+p0)*(d0+3*p0)*TMath::Power(r0,5);
	  d=d0-2/3*3.14*d0/p0*(d0+p0)*(d0+3*p0)*r0*r0;
	  p=p0-TMath::Abs((TMath::Log10(d0)-TMath::Log10(vDensity.at(vPressure.size()-1-dD-1)))/(TMath::Log10(p0)-TMath::Log10(vPressure.at(vPressure.size()-1-dD-1))))*2*3.14*(d0+p0)*(d0+3*p0)*r0*r0/3;
	}
 
	else{
	
	  double mold=m;
	  double rold=r;
	  double dold=d;
	  double pold=p;
	  double rnew,mnew,dnew,pnew;
	  if(d<TMath::Power(10,14.5)*kappa*G/(c*c)) stepsize=10./(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,5.5)*kappa*G/(c*c)) stepsize=1./(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,3.5)*kappa*G/(c*c)) stepsize=0.01/(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,1.5)*kappa*G/(c*c)) stepsize=0.0001/(TMath::Power(kappa,0.5));
	
	  k0=dmdr(rold,mold,pold,vPressure,vEnergyDensity);
	  l0=dpdr(rold,mold,pold,vPressure,vEnergyDensity);

	  if(pold+stepsize*0.5*l0<vPressure.at(0))break;
	  k1= dmdr(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,vPressure,vEnergyDensity);
	  l1= dpdr(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,vPressure,vEnergyDensity);

	  if(pold+stepsize*0.5*l1<vPressure.at(0))break;
	  k2 = dmdr(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,vPressure,vEnergyDensity);
	  l2 = dpdr(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,vPressure,vEnergyDensity);

	  if(pold+stepsize*l2<vPressure.at(0))break;
	  k3 = dmdr(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,vPressure,vEnergyDensity);
	  l3 = dpdr(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,vPressure,vEnergyDensity);

	  rnew=rold+stepsize;
	  mnew=mold+(k0+2.*k1+2.*k2+k3)*stepsize/6.;
	  pnew=pold+(l0+2.*l1+2.*l2+l3)*stepsize/6.;

	  if(pnew<vPressure.at(0))break;
	  int pos=crust(pnew,vPressure,vDensity);

	  double a=(TMath::Log10(vDensity.at(pos))-TMath::Log10(vDensity.at(pos+1)))/(TMath::Log10(vPressure.at(pos))-TMath::Log10(vPressure.at(pos+1)));
	  double b=TMath::Power(10,TMath::Log10(vDensity.at(pos+1))-(TMath::Log10(vDensity.at(pos))-TMath::Log10(vDensity.at(pos+1)))/(TMath::Log10(vPressure.at(pos))-TMath::Log10(vPressure.at(pos+1)))*TMath::Log10(vPressure.at(pos+1)));
	  dnew=b*TMath::Power(pnew,a);
	
	

 
      
	  m=mnew;
	  r=rnew;
	  p=pnew;
	  d=dnew;
	  


	}
	count++;

      }
      //RK4 end
      //cout<<(m*c*c*TMath::Power(kappa,0.5))/(Mo*G)<<endl;
      count=0;
      /*    
      if(dD==0){
	mbefore=m;
	rbefore=r;
      }
      else{
	if(m>mbefore){
	  mbefore=m;
	  rbefore=r;
	}
	if(!boolMax){
	  v_m.push_back((m*c*c*TMath::Power(kappa,0.5))/(Mo*G));
	  v_r.push_back(r*TMath::Power(kappa,0.5)/100000.);
	}
	if(boolMax && m<mbefore){
	  stable_size=dD-1;
	  Maxmass=(mbefore*c*c*TMath::Power(kappa,0.5))/(Mo*G);
	  v_m.push_back((mbefore*c*c*TMath::Power(kappa,0.5))/(Mo*G));
	  v_r.push_back(rbefore*TMath::Power(kappa,0.5)/100000.);
	  v_m.push_back((m*c*c*TMath::Power(kappa,0.5))/(Mo*G));
	  v_r.push_back(r*TMath::Power(kappa,0.5)/100000.);

	  boolMax=false;
	}

      }
      if (Maxmass<Max_mass_limit){
	break;
      }
      */
      
      
      v_m.push_back((m*c*c*TMath::Power(kappa,0.5))/(Mo*G));
      v_r.push_back(r*TMath::Power(kappa,0.5)/100000.);
      count=0;
      if(dD==0){
	Maxmass=m*c*c*TMath::Power(kappa,0.5)/(Mo*G);
      }
      else if(dD>1 && m*c*c*TMath::Power(kappa,0.5)/(Mo*G)>Maxmass){
	Maxmass=0;
	break;
      }
      else if( dD>=1 && m*c*c*TMath::Power(kappa,0.5)/(Mo*G)<Maxmass){
        if(Maxmass<Max_mass_limit) break;
      }
      
    }
    //check the Maximum Mass
    //if(fp.maxMassOrRadius(vDensity,vPressure,vEnergyDensity,0)<Max_mass_limit) continue;

    if(Maxmass<Max_mass_limit) continue;
    

    //cout<<Maxmass<<endl;
    //Go up to real MaxMass

    int counts_to_max=0;
    while(true){
      d=vDensity.at(vDensity.size()-1)+dens_stepsize;
      p=fp.pressure(d);
      v_d.push_back(d*c*c/(G*kappa));
      v_p.push_back(p*c*c*c*c/(G*kappa));
      v_e.push_back(fp.e(d)*c*c/(G*kappa));
	
      vPressure.push_back(p);
      vEnergyDensity.push_back(fp.e(d));
      vDensity.push_back(d);	
      //Finds N by using the politropic equation
      v_h.push_back(TMath::Log(fp.h(d))*(c*c));
      vEnthalpy.push_back(TMath::Log(fp.h(d)));
      vSpeedofSound.push_back(fp.soundSpeedValue(p,d));
      
      if( (v_h.at(v_h.size()-1)<v_h.at(v_h.size()-2))  || (v_p.at(v_p.size()-1)<v_p.at(v_p.size()-2)) || (v_d.at(v_d.size()-1)<v_d.at(v_d.size()-2)) || (v_e.at(v_e.size()-1)<v_e.at(v_e.size()-2))){
	CheckCont=false;
	break;
      }
      BoolSpeed=fp.soundSpeed(p,d);
      up=fp.CheckUpRange(p,d);
      down=fp.CheckDownRange(p,d);
      if(!BoolSpeed)break;
      if(!up)break;
      if(!down)break;

      
      
      double stepsize=100./(TMath::Power(kappa,0.5));
      double r0=TMath::Power(10,-15)/(TMath::Power(kappa,0.5));
      double m0=TMath::Power(10,-15)*G/(c*c*TMath::Power(kappa,0.5));
      double d0=vDensity.at(vDensity.size()-1);//(TMath::Power(10,14.9)+dD*TMath::Power(10,15)*0.1)*kappa*G/(c*c);
      double p0=vPressure.at(vPressure.size()-1);//pressure(d0);
      double maxSoundSpeed=0;

      
      //cout<<d0/(1.66*TMath::Power(10,-24)*kappa*G/(c*c))<<endl;

    
      double r,m,d,p,k0,k1,k2,k3,l0,l1,l2,l3;

      d=d0;
      p=p0;
      m=m0;
  
      //RK4 Begining
      while(d>TMath::Power(10,0.9)*kappa*G/(c*c)){   
	if(count==0){
	  r=r0;
	  m=4*TMath::Pi()*d0*r0*r0*r0/3-TMath::Abs((TMath::Log10(d0)-TMath::Log10(vDensity.size()-1-1))/(TMath::Log10(p0)-TMath::Log10(vPressure.size()-1-1)))*8/15*3.14*3.14*d0/p0*(d0+p0)*(d0+3*p0)*TMath::Power(r0,5);
	  d=d0-2/3*3.14*d0/p0*(d0+p0)*(d0+3*p0)*r0*r0;
	  p=p0-TMath::Abs((TMath::Log10(d0)-TMath::Log10(vDensity.at(vPressure.size()-1-1)))/(TMath::Log10(p0)-TMath::Log10(vPressure.at(vPressure.size()-1-1))))*2*3.14*(d0+p0)*(d0+3*p0)*r0*r0/3;
	}
 
	else{
	
	  double mold=m;
	  double rold=r;
	  double dold=d;
	  double pold=p;
	  double rnew,mnew,dnew,pnew;
	  if(d<TMath::Power(10,14.5)*kappa*G/(c*c)) stepsize=10./(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,5.5)*kappa*G/(c*c)) stepsize=1./(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,3.5)*kappa*G/(c*c)) stepsize=0.01/(TMath::Power(kappa,0.5));
	  if(d<TMath::Power(10,1.5)*kappa*G/(c*c)) stepsize=0.0001/(TMath::Power(kappa,0.5));
	
	  k0=dmdr(rold,mold,pold,vPressure,vEnergyDensity);
	  l0=dpdr(rold,mold,pold,vPressure,vEnergyDensity);

	  if(pold+stepsize*0.5*l0<vPressure.at(0))break;
	  k1= dmdr(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,vPressure,vEnergyDensity);
	  l1= dpdr(rold+stepsize/2.,mold+stepsize*0.5*k0,pold+stepsize*0.5*l0,vPressure,vEnergyDensity);

	  if(pold+stepsize*0.5*l1<vPressure.at(0))break;
	  k2 = dmdr(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,vPressure,vEnergyDensity);
	  l2 = dpdr(rold+stepsize/2.,mold+stepsize*0.5*k1,pold+stepsize*0.5*l1,vPressure,vEnergyDensity);

	  if(pold+stepsize*l2<vPressure.at(0))break;
	  k3 = dmdr(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,vPressure,vEnergyDensity);
	  l3 = dpdr(rold+stepsize,mold+stepsize*k2,pold+stepsize*l2,vPressure,vEnergyDensity);

	  rnew=rold+stepsize;
	  mnew=mold+(k0+2.*k1+2.*k2+k3)*stepsize/6.;
	  pnew=pold+(l0+2.*l1+2.*l2+l3)*stepsize/6.;

	  if(pnew<vPressure.at(0))break;
	  int pos=crust(pnew,vPressure,vDensity);

	  double a=(TMath::Log10(vDensity.at(pos))-TMath::Log10(vDensity.at(pos+1)))/(TMath::Log10(vPressure.at(pos))-TMath::Log10(vPressure.at(pos+1)));
	  double b=TMath::Power(10,TMath::Log10(vDensity.at(pos+1))-(TMath::Log10(vDensity.at(pos))-TMath::Log10(vDensity.at(pos+1)))/(TMath::Log10(vPressure.at(pos))-TMath::Log10(vPressure.at(pos+1)))*TMath::Log10(vPressure.at(pos+1)));
	  dnew=b*TMath::Power(pnew,a);
	
       
	  m=mnew;
	  r=rnew;
	  p=pnew;
	  d=dnew;
	  
	}
	count++;
      }
      //RK4 end
      count=0;
      if((m*c*c*TMath::Power(kappa,0.5))/(Mo*G)>Maxmass){
	v_m.insert(v_m.begin(),(m*c*c*TMath::Power(kappa,0.5))/(Mo*G));
	v_r.insert(v_r.begin(),r*TMath::Power(kappa,0.5)/100000.);
	Maxmass=(m*c*c*TMath::Power(kappa,0.5)/(Mo*G));
      }
      else if((m*c*c*TMath::Power(kappa,0.5))/(Mo*G)<Maxmass){
       	v_d.pop_back();
	v_p.pop_back();
	v_e.pop_back();
	
	vPressure.pop_back();
	vEnergyDensity.pop_back();
	vDensity.pop_back();	
	//Finds N by using the politropic equation
	v_h.pop_back();
	vEnthalpy.pop_back();
	vSpeedofSound.pop_back();
	if(counts_to_max>3){
	  break;
	}
	else{
	  dens_stepsize=dens_stepsize/10.;
	  counts_to_max++;
	}

	     
      }
    }

    //Check the speed of sound
    if(!up) continue;
    if(!down) continue;
    if(!BoolSpeed) continue;
    if(!CheckCont) continue;






   if(v_r.at(0)<13) continue;
    
    
 
    //Write EOS parameters Begining
    string title="";
    title=title+"eosPol"+fileNum;
    eosTable<<title<<" "<<par1<<" "<<par2<<" "<<par3<<" "<<par4<<" "<<fp.rho1*c*c/(kappa*G)<<" "<<fp.rho2*c*c/(kappa*G)<<" "<<surfEOS.at(sEOSnum)<<" "<<endl;
    ofstream output;
    output.open(title,ios::out);
    //Write EOS parameters End

  
    //output Data
    output<<v_d.size()/*-data*/<<endl;//remove data when suurface included
    for(int i=0;i<v_d.size();i++){
      output<<v_e.at(i)<<" "<<v_p.at(i)<<" "<<v_h.at(i)<<" "<<v_d.at(i)/(Mb)<<endl;
    }




    string title2="";
    title2=title2+"eosPol"+fileNum+"table";
    ofstream output2;
    output2.open(title2,ios::out);
    output2<<"Mass Radius SpeedofSound CentralEnergyDensity"<<endl;
    for(int i=0;i<v_d.size()-data;i++){
      output2<<v_m.at(i)<<" "<<v_r.at(i)<<" "<<vSpeedofSound.at(vSpeedofSound.size()-1-i)<<" "<<v_e.at(v_d.size()-1-i)<<endl;
    }

    TGraph* m_r_cur =new TGraph(v_m.size(),&v_r.at(0),&v_m.at(0));
    m_r_cur->SetLineWidth(3.5);
    //m_r->SetLineColor(int(rand->Uniform(0.5,9.5)));
    m_r_cur->SetLineColor(fileNum+1);
    m_r_cur->SetMarkerColor(fileNum+1);
    m_r_cur->SetMarkerStyle(20);
    m_r_cur->SetMarkerSize(0.5);
    m_r_cur->Draw("L SAME");
    canvas->SaveAs("Mass_Radius.png");

    fileNum++;
    //reset the vectors
    output.close();
    output2.close();
    vEnergyDensity.clear();
    vEnthalpy.clear();
    vDensity.clear();
    vPressure.clear();
  }
  //Random EOS Generation Begins

 

}
