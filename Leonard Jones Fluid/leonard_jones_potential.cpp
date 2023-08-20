#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include <fstream>
//#include <sciplot/sciplot.hpp>

#define myrand() (rand()/double(RAND_MAX))

//using namespace sciplot;
using namespace std;

double energy(double * xx, double *yy,double &vir,int N,double L) // Engergy calculation
{	double r2,r6,e=0.0,dx,dy,hL=L/2.0;
	vir =0.0;
	for(int i=0; i<N; i++)
	{
		for(int j=i+1;j<N;j++)
		{
			dx = xx[i]-xx[j];
			dy = yy[i]-yy[j];


			if (dx>hL)       dx-=L; //boundary conditions
			else if (dx<-hL) dx+=L;
			if (dy>hL)       dy-=L;
			else if (dy<-hL) dy+=L;

			r2 = dx*dx + dy*dy; // r is the distance
			r6 = 1.0/(r2*r2*r2); // calculating r^6 
			e += 4*(r6*r6 - r6);
			vir += 24*(r6*r6-0.5*r6);

		}
	}
	return e;
}


void init(double * xx, double *yy, double L,int N) //to initialize the configuration
{	int c = 0,n=sqrt(N);
	cout<< N<<" "<<n<<endl;
	cout<< L<<endl;
	double d= L/n; //lattice constant 
	cout<< d<<endl;
	//exit(0);
	for(int i=0; i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			xx[c]=j*d;
			yy[c]=i*d;
			c++;		
		}
	}
}


double metropolis(double * xx,double *yy,double V,double rho,int cycle,int N,double L,double T,ofstream& fp,double *p_c)
{
	int pr;
	double x_old, y_old, E_old, E_new,dr=0.1,dx,dy;
	double racc,facc=0.0;
	double vir_new=0.0,vir_old=0.0,vir_sum=0.0,p,p_sum=0.0;
	E_old = energy(xx,yy,vir_old,N,L);
	
	ofstream simf("LJ_sim.txt");

	//cout<< "OKAY1"<<endl;
	int count=0;
	for(int c=0;c<cycle; c++)
	{
		
		for(int id=0;id<N;id++)
		{
		pr = myrand()*N; //randomly choose a particle by its index
		
		
		// Record the previous positions
		x_old=xx[pr];
		y_old=yy[pr];

		dx = dr*(2*myrand()-1);
		dy = dr*(2*myrand()-1);

		xx[pr] += dx;
		yy[pr] += dy;

		//Apply periodic boundary conditions
		if(xx[pr]<0.0) xx[pr]+=L;
		if(xx[pr]>L) xx[pr]-=L;
		if(yy[pr]<0.0) yy[pr]+=L;
		if(yy[pr]>L) yy[pr]-=L;


		E_new = energy(xx,yy,vir_new,N,L);

		if(myrand()<exp(-T*(E_new-E_old)))
		{
			E_old=E_new;
			vir_old = vir_new;
			facc +=1;
		}
		else
		{
			xx[pr]=x_old;
			yy[pr]=y_old;
		}

		//vir_sum += vir_old;
		simf << xx[id]<<" "<<yy[id]<<"\n";
		}
		// Adjusting maximum displacement according to the system parameters
		racc = facc/N;
		if (racc >0.5)
		{dx *= 1.05;dy *= 1.05;}
		else {dx *= 0.95;dy *= 0.95;}

		simf <<"\n"<<"\n";
		
		p = rho*(1/T)+vir_old/V;
		p_c[c] = p;
	// Cycle vs Energy plot
	fp << c<<" "<<(E_old/N)<<" "<<p<<endl;
	cout << c << " "<<(E_old/N)<<" "<<p<<endl;
	
	//Plot plot;
	//plot.drawCurve(c,p);
	
	
	if(c>5000 && c%100==0)
	{
		p_sum += p;count++;
	} 


	}
	simf.close();

	return (p_sum/(cycle-5000-count));
	
	


	// Experimental code for ploting using Sciplot
	//vector<double> x_pos(xx,xx+sizeof(xx)/sizeof(double));

	//vector<double> y_pos(yy,yy+sizeof(yy)/sizeof(double));
	//vector<double> y_pos(std::begin(yy),std::end(yy)); 
	//for(int i:x_pos){
	//cout<<i<<endl;}
	//Plot plot;
	//plot.palette("set2");
	//plot.drawCurve(x_pos,y_pos);
	//plot.show();

}

double correlation(double *p,int N, int t)
{
	//C(t) = (<s(T)s(T+t)> - <s(T)>)/(<s^2(T)> - <s(T)>^2);
	double p_sum=0,p2_sum=0,pt_sum=0,p_avg,cor;
	for(int i=0; i<N; i++)
	{
		p_sum += p[i];
		p2_sum += p[i]*p[i];
		if ((i+t) > N) continue;
		pt_sum += p[i]*p[i+t];

	}
	p_avg = p_sum/N;
	cor = (pt_sum/(N-t) -p_avg*p_avg)/(p2_sum/N - p_avg*p_avg);
	
	return cor;
}


int main()
{
 
 int cycle=10000;
 int N=10*10;  // NUmber of particles = 100
 double T,V,L,rho,p_avg;
 double xx[N],yy[N];
 double p_c[cycle];//array to store pressure so to calculate correlation afterwards


 cout << "Enter the temperature"<<endl;
 //cin >> T;
 T = 2;
 T = 1/T;// reciprocal of Temperature for efficiency
 rho = 0.4;
 V = N/rho; //Volume = N/rho 
 L = sqrt(V);
 
 
 ofstream fp("Energy.txt");

 init(xx,yy,L,N);
 cout<<"OKAY"<<endl;
 p_avg = metropolis(xx,yy,V,rho,cycle,N,L,T,fp,p_c);

 fp.close();
 
 cout<<p_avg<<" "<<rho<<endl;



 FILE* gnuplot;
 gnuplot = popen("gnuplot -p", "w");
 if(gnuplot != NULL){fprintf(gnuplot, "p 'Energy.txt' usi 1:3 w l\n");}
 fprintf(gnuplot, "p 'Energy.txt' usi 1:2 w l \n");
 //fprintf(gnuplot, "p 'Energy.txt' usi 1:3\n");
 
 system("gnuplot LJ_simulation.gnu");
 ofstream fcor("correlation.txt");
 
 double cor;
 for(int i=1; i<100; i++)
 {
 fcor<<i<<" "<<correlation(p_c,N,i)<<endl;
 }

 fcor.close();
}

