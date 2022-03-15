#define n0 0.16*TMath::Power(10,39)//1/cm^3
#define mN 1.675*TMath::Power(10,-24)//gr
#define G 6.6732*TMath::Power(10,-8)
#define c 2.9979*TMath::Power(10,10)
#define Mo 1.987*TMath::Power(10,33)
#define eo  TMath::Power(10,15)
#define Mb  1.66*TMath::Power(10,-24)
#define kappa c*c/(G*eo)
#define hbar 1.054571628*TMath::Power(10,-34)





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
  


double f(double e,double a1,double a2,double a3,double a4,double a5,double a6){
  
  double speed;
  //if(e/(mN*n0)>1.5){
    speed=a1*TMath::Exp(-0.5*TMath::Power(e/(mN*n0)-a2,2)/(a3*a3))+a6+(1./3.-a6)/(1.+TMath::Exp(-a5*((e/(mN*n0)-a4))));
  //}
  //else{
    //speed=(1-0.52)*(hbar*hbar)*TMath::Power(3.*TMath::Power(TMath::Pi(),2)*e/(mN),2./3.)/(0.92*3.*mN*mN*c*c);
    //speed=(1-0.52)*TMath::Power(3.*TMath::Power(TMath::Pi(),2)*e/(mN),2./3.)/(0.92*3.*2.2676*TMath::Power(10,27));
  //}


  if (speed<0){
    return 0;
  }
  else{
    return speed;
  }
}

double pressure(double e,double e0,double p0,double a1,double a2,double a3,double a4,double a5,double a6){
  double dx=(e-e0)/1000.;
  double integ=0;
  double ener=e0;
  integ=f(e0,a1,a2,a3,a4,a5,a6)+f(e,a1,a2,a3,a4,a5,a6);

    
  for(int i=1;i<1000-1;i++){
    ener+=dx;
    if(i%2==0){
      integ+=2.*f(ener,a1,a2,a3,a4,a5,a6);
    }
    else{
      integ+=4.*f(ener,a1,a2,a3,a4,a5,a6);
    }
  }

  return p0+integ*dx*c*c/3.;
}

double g(double e,double e0, double p0,double a1,double a2,double a3,double a4,double a5,double a6){
  
  return 1./(e+pressure(e,e0,p0,a1,a2,a3,a4,a5,a6)/(c*c));
}


double enthalpy(double e,double e0,double p0,double h0,double a1,double a2,double a3,double a4,double a5,double a6){
  double dx=(e-e0)/100.;
  double integ=0;
  double ener=e0;
  integ=g(e0,e0,p0,a1,a2,a3,a4,a5,a6)+g(e,e0,p0,a1,a2,a3,a4,a5,a6);

    
  for(int i=1;i<100-1;i++){
    ener+=dx;
    if(i%2==0){
      integ+=2.*g(ener,e0,p0,a1,a2,a3,a4,a5,a6);
    }
    else{
      integ+=4.*g(ener,e0,p0,a1,a2,a3,a4,a5,a6);
    }
  }
  return h0+integ*(pressure(e,e0,p0,a1,a2,a3,a4,a5,a6)-pressure(e0,e0,p0,a1,a2,a3,a4,a5,a6))/(3*100.);
}




void greif(){

  //asks for EOS num
  int num;
  cout<<"Give the number of EOS"<<endl;
  cin>>num;

  int interv=40;
 
  double a1,a2,a3,a4,a5,a6;

  ofstream eosTable;
  eosTable.open("EosTable",ios::out);
  eosTable<<"EOS  a1  a2  a3  a4  a5  a6  central_energy cEFT"<<endl;
  
  TRandom3* rand =new TRandom3();
  rand->SetSeed(0);
  

  double vecA1[12]={0.533484,0.806516,0.358644,1.31031,1.39068,1.32732,1.46422,0.931557,0.686456,0.824155,1.16696,1.1814};
  double vecA2[12]={6.02743,4.27762,2.18641,4.07643,3.35555,3.12917,6.12736,2.99786,5.34602,3.62393,4.26464,4.87925};
  double vecA3[12]={4.70183,1.79287,2.9905,2.11916,1.89917,2.35747,5.53277,0.846201,2.20362,0.853554,4.03632,3.136};
  double vecA4[12]={2.46758,8.78623,2.19397,9.14079,6.94482,4.94596,12.4406,2.95033,9.49444,21.1719,6.80565,6.82788};
  double vecA5[12]={0.829607,0.75933,0.961889,0.497626,0.860675,0.951592,0.495523,0.417567,0.54512,0.620857,0.787862,0.384679};

  int fileNum=0;


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


  
  //Loop for EOS
  while(fileNum<num){
    double e0=1.*TMath::Power(10,15); /* kappa*G/(c*c) *///central energy density

    vector<double> v_p,v_d,v_e,v_h,v_r,v_m;
    vector<double> vEnergyDensity,vEnthalpy,vDensity,vPressure; //non dimensional
    
    //choose surface EOS
    vector<string> surfEOS;
    surfEOS.push_back("eosNV");
    surfEOS.push_back("eosFPS");
    surfEOS.push_back("eosBPS");
    int sEOSnum=2;//0 = NV   1 = FPS 2 = BBB2
    //Maximum Mass
    double Max_mass_limit=1.97; //solar masses




    ifstream myinputfile;
    myinputfile.open(surfEOS.at(sEOSnum),ios::in);
    string line;
    string read1,read2,read3,read4,read5;
    int data;
    int count=0;
    while( getline(myinputfile,line) ){
	stringstream liness = stringstream(line);
	getline(liness,read1,' ');
	getline(liness,read2,' ');
	getline(liness,read3,' ');
	getline(liness,read4,' ');
	if(sEOSnum==0){getline(liness,read5,' ');}
	if(stod(read4)>0.5*n0) break;
	v_e.push_back(stod(read1));
	v_p.push_back(stod(read2));
	v_h.push_back(stod(read3));
	v_d.push_back(stod(read4)*Mb);
	vPressure.push_back(stod(read2)*kappa*G/(c*c*c*c));
	vEnergyDensity.push_back(stod(read1)*kappa*G/(c*c));
	vEnthalpy.push_back(stod(read3)/(c*c));
	vDensity.push_back(stod(read4)*Mb*kappa*G/(c*c));
    }
    data=v_d.size();

     
    //cEFT band
    string cEFT;
    if(rand->Uniform(0,1)>0.5){
    //stiff
    cEFT="eoscEFTstiff";
    }
    else{
    //soft
    
    cEFT="eoscEFTsoft";
    }
    
    ifstream myinputfile1;
    myinputfile1.open(cEFT,ios::in);


  count=0;
    while( getline(myinputfile1,line) ){
	stringstream liness = stringstream(line);
	getline(liness,read1,' ');
	getline(liness,read2,' ');
	getline(liness,read3,' ');
	getline(liness,read4,' ');

	v_e.push_back(stod(read1));
	v_p.push_back(stod(read2));
	v_h.push_back(stod(read3));
	v_d.push_back(stod(read4)*Mb);
	vPressure.push_back(stod(read2)*kappa*G/(c*c*c*c));
	vEnergyDensity.push_back(stod(read1)*kappa*G/(c*c));
	vEnthalpy.push_back(stod(read3)/(c*c));
	vDensity.push_back(stod(read4)*Mb*kappa*G/(c*c));
    }    
    data=v_d.size();
     
          
    double e=v_e.at(v_e.size()-1);
    double de=(e0-e)/double(interv);      
    bool BoolSpeed=true;

    a1=rand->Uniform(0.1,1.5);//vecA1[fileNum];//0.533484;//
    a2=rand->Uniform(1.5,12);//vecA2[fileNum];// 6.02743;//
    a3=rand->Uniform(0.05*a2,2*a2);//vecA3[fileNum];//4.70183;//
    a4=rand->Uniform(1.5,37);//vecA4[fileNum];//2.46758;//
    a5=rand->Uniform(0.1,1);//vecA5[fileNum];//0.829607;//

    double u =TMath::Sqrt((v_p.at(v_p.size()-1)-v_p.at(v_p.size()-2))/(v_e.at(v_e.size()-1)-v_e.at(v_e.size()-2)))/(c);

    a6=/*-0.0395129;*/(u*u-a1*TMath::Exp(-0.5*TMath::Power(v_e.at(v_e.size()-1)/(mN*n0)-a2,2)/(a3*a3))-1./(3.*(1.+TMath::Exp(-a5*((v_e.at(v_e.size()-1)/(mN*n0)-a4))))))/(1.-1./(1.+TMath::Exp(-a5*((v_e.at(v_e.size()-1)/(mN*n0)-a4)))));
    //a6=/*-0.0395129;*/(0.05-a1*TMath::Exp(-0.5*TMath::Power(TMath::Power(10,14)/(mN*n0)-a2,2)/(a3*a3))-1./(3.*(1.+TMath::Exp(-a5*(3.*TMath::Power(10,14)/(mN*n0)-a4)))))/(1.-1./(1.+TMath::Exp(-a5*((TMath::Power(10,14)/(mN*n0)-a4)))));
    
    //cout<<u<<endl;
    if(f(TMath::Power(10,16),a1,a2,a3,a4,a5,a6)>(1./3.+0.001) || f(TMath::Power(10,16),a1,a2,a3,a4,a5,a6)<(1./3.-0.001)) continue;

    for(int i=0;i<interv;i++){

      e+=de;
      double pr=pressure(e,v_e.at(v_e.size()-1),v_p.at(v_p.size()-1),a1,a2,a3,a4,a5,a6);
      double h=enthalpy(e,v_e.at(v_e.size()-1),v_p.at(v_p.size()-1),v_h.at(v_h.size()-1),a1,a2,a3,a4,a5,a6);
      double n= (e+pr/(c*c))*exp(-h/(c*c))/(Mb);
      //cout<<e<<" "<<pr<<" "<<f(e,a1,a2,a3,a4,a5,a6)<<" "<<TMath::Sqrt((pr-v_p.at(v_p.size()-1))/(e-v_e.at(v_e.size()-1)))/(c)<<endl;
      v_p.push_back(pr);
      v_e.push_back(e);
      v_d.push_back(n*Mb);
      v_h.push_back(h);
      vPressure.push_back(pr*kappa*G/(c*c*c*c));
      vEnergyDensity.push_back(e*kappa*G/(c*c));
      vDensity.push_back(n*Mb*kappa*G/(c*c));
      vEnthalpy.push_back(h/(c*c));
      
      if(f(e,a1,a2,a3,a4,a5,a6)>1 || f(e,a1,a2,a3,a4,a5,a6)<=0){
	BoolSpeed=false;
	break;
      }
      if(n<1.5*n0 && f(e,a1,a2,a3,a4,a5,a6)>0.163){
	BoolSpeed=false;
	break;
      }

    }
    


    double dume=e0;

    while(1){
      if(f(dume,a1,a2,a3,a4,a5,a6)>1 || f(dume,a1,a2,a3,a4,a5,a6)<=0){
	BoolSpeed=false;
	break;
      } 
    if(dume>1.*TMath::Power(10,16)) break;
    dume+=de;
       // cout<<dume<<endl;
    }


    if(!BoolSpeed) continue;
    //if(maxMassOrRadius(vDensity,vPressure,vEnergyDensity,0)<1.97) continue;


    double Maxmass=0;
    for(int dD=0;dD<interv-5;dD++){

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
      v_m.push_back((m*c*c*TMath::Power(kappa,0.5))/(Mo*G));
      v_r.push_back(r*TMath::Power(kappa,0.5)/100000.);
      count=0;
      if(dD==0){
	Maxmass=m*c*c*TMath::Power(kappa,0.5)/(Mo*G);
      }
      else if(dD>1 && m*c*c*TMath::Power(kappa,0.5)/(Mo*G)>Maxmass){
	Maxmass=m*c*c*TMath::Power(kappa,0.5)/(Mo*G);
	//break;
      }
      //else if( dD>=1 && m*c*c*TMath::Power(kappa,0.5)/(Mo*G)<Maxmass){
        //if(Maxmass<Max_mass_limit) break;
      //}
      cout<<m*c*c*TMath::Power(kappa,0.5)/(Mo*G)<<endl;
    }
    //cout<<Maxmass<<endl;
    //check the Maximum Mass
    //if(fp.maxMassOrRadius(vDensity,vPressure,vEnergyDensity,0)<Max_mass_limit) continue;

    //if(Maxmass<1.97) continue;




    int counts_to_max=0;
    while(true){
      e=v_e.at(v_e.size()-1)+de;
      double pr=pressure(e,v_e.at(v_e.size()-1),v_p.at(v_p.size()-1),a1,a2,a3,a4,a5,a6);
      double h=enthalpy(e,v_e.at(v_e.size()-1),v_p.at(v_p.size()-1),v_h.at(v_h.size()-1),a1,a2,a3,a4,a5,a6);
      double n= (e+pr/(c*c))*exp(-h/(c*c))/(Mb);
      //cout<<e<<" "<<pr<<" "<<f(e,a1,a2,a3,a4,a5,a6)<<" "<<TMath::Sqrt((pr-v_p.at(v_p.size()-1))/(e-v_e.at(v_e.size()-1)))/(c)<<endl;
      v_p.push_back(pr);
      v_e.push_back(e);
      v_d.push_back(n*Mb);
      v_h.push_back(h);
      vPressure.push_back(pr*kappa*G/(c*c*c*c));
      vEnergyDensity.push_back(e*kappa*G/(c*c));
      vDensity.push_back(n*Mb*kappa*G/(c*c));
      vEnthalpy.push_back(h/(c*c));
      
      if(f(e,a1,a2,a3,a4,a5,a6)>1 || f(e,a1,a2,a3,a4,a5,a6)<=0){
	BoolSpeed=false;
	break;
      }
      if((e/(mN*n0)<1.5) && f(e,a1,a2,a3,a4,a5,a6)>0.163){
	BoolSpeed=false;
	break;
      }


    
    if(!BoolSpeed) continue;
      
      
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
      cout<<m*c*c*TMath::Power(kappa,0.5)/(Mo*G)<<endl;

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
	if(counts_to_max>3){
	  break;
	}
	else{
	  de=de/10.;
	  counts_to_max++;
	}
    //cout<<Maxmass<<endl;
	     
      }

    }

    if(Maxmass<Max_mass_limit) continue;
    if(!BoolSpeed) continue;








    //Write EOS parameters Begining
    string title="";
    title=title+"eosGreif"+fileNum;
    eosTable<<title<<" "<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<v_e.at(v_e.size()-1)<<" "<<cEFT<<endl;
    
    ofstream output;
    output.open(title,ios::out);
    //Write EOS parameters End
    //output Data
    output<<v_d.size()<<endl;
    for(int i=0;i<v_d.size();i++){
      output<<v_e.at(i)<<" "<<v_p.at(i)<<" "<<v_h.at(i)<<" "<<v_d.at(i)/(Mb)<<endl;
    }


    string title2="";
    title2=title2+"eosGreif"+fileNum+"table";
    ofstream output2;
    output2.open(title2,ios::out);
    output2<<"Mass Radius CentralEnergyDensity"<<endl;
    for(int i=0;i<v_m.size();i++){
      output2<<v_m.at(i)<<" "<<v_r.at(i)<<" "<<v_e.at(v_e.size()-1-i)<<endl;
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

    output.close();
    v_e.clear();
    v_p.clear();
    v_d.clear();
    v_h.clear();
    fileNum++;
  }
}
