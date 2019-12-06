#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

using namespace std;

//Global Variables
int i,j,nh=128,nv=128,n=nv+2;
double dx=1.0/n;
double dy=1.0/n;
double dt = 0.00001;
float Re = 100.0;
double div = 0.0;
double limit = 1e-8;
double limit_temp=1e-6;
int inner_loop;
double t=0.0;
double error_U =0.0;
double error_V =0.0;
double error_P=0.0;
double error_T = 0.0;
double error_t=0.0;
double s=0.0;
double s_U=0.0;
double s_V=0.0;
double s_T=0.0;
double S_s =0.0;
double t_t =0.0;
float Pr=0.71;
float Ra= 1e6;
double nu_wall_west=0.0;
double nu_wall_east=0.0;
float w=1.0;
int outer_loop=0;


//declearing all subroutines
void grid_generation();
void grid_output();
void equation_3b();
void equation_3a();
void Velocity_BC();
void equation_1a();
void equation_1b();
void Y_momentum_solver();
void X_momentum_solver();
void intialisation();
void grid_output();
void MAC();
void correct_pressure();
void pressure_solver();
void pressure_coffecients();
void swap();
void store_old_velocities();
void velocity_divergence();
void residual_print();
void P_tilda_store();
void swap_p_tilda();
void U_plot();
void V_plot();
void CT();
void DT();
void temp();
void temp_boundary_conditions();
void T_plot();
void vector_plot();
void Nusselt_plot();




struct data{
		double P_x,P_y,U_x,U_y,V_x,V_y;
        double U,V,P,U_old,V_old,P_old;
        double ue2,uw2,du2dx,un,vn,us,vs,duvdy,Cu;
        double ue2_old,uw2_old,du2dx_old,un_old,vn_old,us_old,vs_old,duvdy_old,Cu_old;
        double d2udx2,d2udy2,Du;
        double d2udx2_old,d2udy2_old,Du_old;
        double vn2,vs2,dv2dy,ue,ve,uw,vw,duvdx,Cv;
        double vn2_old,vs2_old,dv2dy_old,ue_old,ve_old,uw_old,vw_old,duvdx_old,Cv_old;
        double d2vdx2,d2vdy2,Dv;
        double d2vdx2_old,d2vdy2_old,Dv_old;
        double U_tilda,V_tilda,P_tilda;
        double rho_temp,rho;
        double domain_x,domain_y;
        double AP,AS,AW,AN,AE,BP,P_tilda_new,P_tilda_old;
        double DIV;
        double P_new;
        double T,T_new,T_old;
        double Ct,uedtedx,uwdtwdx,vndtndy,vsdtsdy;
        double Dt,d2tedx2,d2twdx2,d2tndy2,d2tsdy2;
        double Nu;

        }node[200][200];

int main()
{
	clock_t begin, end;
    double cpu_time_used;
    begin = clock();

    grid_generation();
    intialisation();
    grid_output();
    MAC();
    end = clock();
	cpu_time_used = ((double) (end - begin))/ CLOCKS_PER_SEC;
	cout<<" \n"<<cpu_time_used;
}
void MAC(){
		do{ outer_loop++;
			Velocity_BC();
          	store_old_velocities();
			equation_1a();//pedicting u velocity
			equation_1b();//predicting v velocity
			Velocity_BC();
			velocity_divergence();
			//pressure_coffecients();
					inner_loop=0;
				do{
					inner_loop++;
				//	pressure_solver();
					swap_p_tilda();
					equation_3a();//correcting u velocity
					equation_3b();//correcting v velocity
					velocity_divergence();
					error_P = sqrt(s/pow(n,2));

				//	cout<<inner_loop<<"\t"<<error_P<<"\n";
				}while(error_P>limit);
					correct_pressure();
					swap();

					t=t+dt;

				do{

					temp();
					error_t = sqrt(t_t/pow(n,2));
					cout<<"temp error"<<error_t<<"\n";
				}while(error_t>limit_temp);
					error_U = sqrt(s_U/pow(n,2));
					error_V = sqrt(s_V/pow(n,2));
					error_T = sqrt(s_T/pow(n,2));
					residual_print();
					cout<<"uerror = "<<error_U<<"\t"<<"verror = "<<error_V<<" pressure "<<S_s<<" temperature "<<error_T<<"\n";
					cout<<"time = "<<t<<"\n";
					cout<<"inner loop time = "<<inner_loop<<"\n";


                         if(outer_loop == 2500 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}

                        else if(outer_loop == 5000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 7500 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 10000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 15000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 20000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();
							}
			else if(outer_loop == 25000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}	
			else if(outer_loop == 30000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 32000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 34000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 36000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 38000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 40000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 42000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 44000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 46000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}
			else if(outer_loop == 48000 ){
								U_plot();
								V_plot();
								T_plot();
								vector_plot();
								Nusselt_plot();

							}


			
				}while(t<=0.5);
				U_plot();
				V_plot();
				T_plot();
				vector_plot();
				Nusselt_plot();
}
void intialisation(){
    for(j=0;j<=n+1;j++){
        for(i=0;i<=n+1;i++){
            node[i][j].U = 0.0;//U initialised
			node[i][j].U_old = 0.0;//U initialised
            node[i][j].V = 0.0;// V initialised
			node[i][j].V_old = 0.0;// V initialised
            node[i][j].P = 1.0;// P initialised
			node[i][j].P_old = 0.0;//pressure correction
			node[i][j].P_tilda_old = 0.0;
			node[i][j].T = 0.0;
			node[i][j].T_new =1.0;
			node[i][j].T_old =1.0;
        }
    }
}
void velocity_divergence(){
	P_tilda_store();

 double a_prime = -2.0*dt*(1.0/pow(dx,2) + 1.0/pow(dy,2) );
                s = 0.0;
    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
                node[i][j].P_tilda_new = w*((node[i][j].U_tilda - node[i-1][j].U_tilda)/dx + (node[i][j].V_tilda - node[i][j-1].V_tilda)/dy)/a_prime;
				s = s+ pow((node[i][j].P_tilda_new-node[i][j].P_tilda_old), 2);
                node[i][j].P_tilda = node[i][j].P_tilda_new;

        }
    }
}
void swap_p_tilda(){
	for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
			node[i][j].P_tilda = node[i][j].P_tilda_new;
			}
		}
}
void X_momentum_solver(){
    //X momentum solved in ucell

    // convective term
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].ue2 =pow(node[i+1][j].U + node[i][j].U,2);
                node[i][j].uw2 = pow(node[i-1][j].U + node[i][j].U,2);
                node[i][j].du2dx = (node[i][j].ue2 - node[i][j].uw2)/(4*dx);

                node[i][j].un = (node[i][j+1].U + node[i][j].U);
                node[i][j].us = (node[i][j-1].U + node[i][j].U);
                node[i][j].vn = (node[i+1][j].V + node[i][j].V);
                node[i][j].vs = (node[i][j-1].V + node[i+1][j-1].V);
                node[i][j].duvdy = (node[i][j].un*node[i][j].vn - node[i][j].us*node[i][j].vs)/(4*dy);

                node[i][j].Cu = node[i][j].du2dx + node[i][j].duvdy;

        }
    }
    //Diffusive term
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].d2udx2 = (node[i+1][j].U - 2.0*node[i][j].U + node[i-1][j].U)/(dx*dx);
                node[i][j].d2udy2 = (node[i][j+1].U - 2.0*node[i][j].U + node[i][j-1].U)/(dy*dy);

                node[i][j].Du = node[i][j].d2udx2 + node[i][j].d2udy2;
        }
    }

}

void Y_momentum_solver(){
    //Y momentum solved in V cell

    // convective term
    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
            node[i][j].vn2 = pow(node[i][j+1].V + node[i][j].V,2);
            node[i][j].vs2 = pow(node[i][j-1].V + node[i][j].V,2);
            node[i][j].dv2dy = (node[i][j].vn2 -node[i][j].vs2)/(4*dy);

            node[i][j].ue = (node[i][j].U + node[i][j+1].U);
            node[i][j].uw = (node[i-1][j].U + node[i-1][j+1].U);
            node[i][j].ve = (node[i][j].V + node[i+1][j].V);
            node[i][j].vw = (node[i-1][j].V + node[i][j].V);
            node[i][j].duvdx = (node[i][j].ue*node[i][j].ve - node[i][j].uw*node[i][j].vw)/(4*dx);

            node[i][j].Cv = node[i][j].dv2dy + node[i][j].duvdx;
        }
    }
    //diffusive term
    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
                node[i][j].d2vdx2 = (node[i+1][j].V - 2.0*node[i][j].V + node[i-1][j].V)/(dx*dx);
                node[i][j].d2vdy2 = (node[i][j+1].V - 2.0*node[i][j].V + node[i][j-1].V)/(dy*dy);

                node[i][j].Dv = node[i][j].d2vdx2 + node[i][j].d2vdy2;
        }
    }


}

//predictor step

//U momentum in U cell
void equation_1a(){

	X_momentum_solver();
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].U_tilda = node[i][j].U - dt*node[i][j].Cu + Pr*dt*node[i][j].Du - (node[i+1][j].P - node[i][j].P)*(dt/dx);
        }
    }

}
//V momentum in V cell

void equation_1b(){

	Y_momentum_solver();
    for(j=1;j<=n-1;j++){
            for(i=1;i<=n;i++){
                    node[i][j].V_tilda = node[i][j].V - dt*node[i][j].Cv + Pr*dt*node[i][j].Dv - (node[i][j+1].P - node[i][j].P)*(dt/dy)+ 0.5*dt*Ra*Pr*(node[i][j].T+node[i][j+1].T);
            }
    }

}
//boundary conditions of velocity
void Velocity_BC(){


    //west wall
    for(j=0;j<=n;j++){
        node[0][j].U = 0.0;//no slip
		node[0][j].U_old = 0.0;//no slip
        node[0][j].U_tilda = 0.0;
        node[0][j].V = -node[1][j].V;
		node[0][j].V_old = -node[1][j].V_old;
        node[0][j].V_tilda = -node[1][j].V_tilda;//no slip
    }
    //east wall
     for(j=0;j<=n;j++){
        node[n][j].U = 0.0;//no slip
		node[n][j].U_old = 0.0;//no slip
        node[n][j].U_tilda = 0.0;//no slip
        node[n+1][j].V = -node[n][j].V;//no slip
		node[n+1][j].V_old = -node[n][j].V_old;//no slip
        node[n+1][j].V_tilda = -node[n][j].V_tilda;
     }
     //top wall
     for(i=0;i<=n;i++){
        node[i][n+1].U =  - node[i][n].U;//no slip
		node[i][n+1].U_old = - node[i][n].U_old;//no slip
        node[i][n+1].U_tilda =- node[i][n].U_tilda;
        node[i][n].V = 0.0;//no slip
		node[i][n].V_old = 0.0;//no slip
        node[i][n].V_tilda = 0.0;//no slip
     }
     //bottom wall
     for(i=0;i<=n;i++){
        node[i][0].U = - node[i][1].U;//no slip
		node[i][0].U_old = - node[i][1].U;//no slip
        node[i][0].U_tilda = - node[i][1].U_tilda;//no slip
        node[i][0].V = 0.0;//no slip
		node[i][0].V_old = 0.0;//no slip
        node[i][0].V_tilda = 0.0;//no slip
     }

}
// u correction in u cell
void equation_3a(){

    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].U_tilda = node[i][j].U_tilda - (dt/dx)*(node[i+1][j].P_tilda - node[i][j].P_tilda);
        }
    }

}
// v correction in v cell
void equation_3b(){

    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
                node[i][j].V_tilda = node[i][j].V_tilda - (dt/dy)*(node[i][j+1].P_tilda - node[i][j].P_tilda);
        }
    }

}
void correct_pressure(){
   swap_p_tilda();
   S_s=0.0;
	for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
			node[i][j].P_new = node[i][j].P + node[i][j].P_tilda;
			S_s = S_s + (pow(node[i][j].P_new - node[i][j].P_old,2));
			node[i][j].P = node[i][j].P_new;
		}
	}
}


void pressure_coffecients(){
    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
            node[i][j].AS = 0.0;
            node[i][j].AN = 0.0;
            node[i][j].AE = 0.0;
            node[i][j].AW = 0.0;
        }
    }
    for(j=2;j<=n;j++){
        for(i=1;i<=n;i++){
            node[i][j].AS = 1.0/(dy*dy);
        }
    }

    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
            node[i][j].AN = 1.0/(dy*dy);
        }
    }

    for(j=1;j<=n;j++){
        for(i=2;i<=n;i++){
            node[i][j].AW = 1.0/(dx*dx);
        }
    }

    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
            node[i][j].AE = 1.0/(dy*dy);
        }
    }
    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
            node[i][j].AP = -(node[i][j].AW + node[i][j].AS + node[i][j].AN + node[i][j].AE);
        }
    }
}

void pressure_solver(){


	s=0;
    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
        	 node[i][j].BP = (1.0/dt)*(node[i][j].U_tilda - node[i-1][j].U_tilda)/dx + (node[i][j].V_tilda - node[i][j-1].V_tilda)/dy;
            node[i][j].rho = node[i][j].BP - (node[i][j].AW*node[i-1][j].P_tilda + node[i][j].AE*node[i+1][j].P_tilda + node[i][j].AN*node[i][j+1].P_tilda + node[i][j].AS*node[i][j-1].P_tilda + node[i][j].AP*node[i][j].P_tilda);
            node[i][j].P_tilda_new = node[i][j].P_tilda + node[i][j].rho/node[i][j].AP;
                s = s+ pow((node[i][j].P_tilda_new-node[i][j].P_tilda), 2);
                node[i][j].P_tilda = node[i][j].P_tilda_new;
        }
    }
}
void swap(){
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].U = node[i][j].U_tilda;
        }
    }
	    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
                node[i][j].V = node[i][j].V_tilda ;
        }
    }
}
void store_old_velocities(){
	for(j=1;j<=n;j++){
		for(i=1;i<=n-1;i++){
		node[i][j].U_old = node[i][j].U;

		}
	}

	for(j=1;j<=n-1;j++){
		for(i=1;i<=n;i++){
		node[i][j].V_old = node[i][j].V;
		}
	}
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){
			node[i][j].P_old = node[i][j].P;
			node[i][j].T_old = node[i][j].T;
			}
		}
}

void grid_generation(){

	for(j=0;j<=n;j++){
		for(i=0;i<=n;i++){
			node[i][j].domain_x = i*dx;
			node[i][j].domain_y = j*dy;
		}
	}

	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){
			node[i][j].P_x = dx/2.0 + (i-1)*dx;
			node[i][j].P_y = dy/2.0 + (j-1)*dy;
			}
		}
	for(j=1;j<=n;j++){
		for(i=0;i<=n;i++){
			node[i][j].U_x = i*dx;
			node[i][j].U_y = dy/2.0 + (j-1)*dy;
			}
		}
	for(j=0;j<=n;j++){
		for(i=1;i<=n;i++){
			node[i][j].V_x = dx/2.0 +(i-1)*dx;
			node[i][j].V_y = j*dy;
			}
		}
}

void grid_output(){
	ofstream f1;
	f1.open("grid.dat");

	f1<<"zone "<<"i= "<<n+1<<" j= "<<n+1<<"\n";
		for(j=0;j<=n;j++){
			for(i=0;i<=n;i++){
				f1<<node[i][j].domain_x<<"\t"<<node[i][j].domain_y<<"\t"<<"1"<<"\n";
				}
			}
	f1<<"zone "<<"i= "<<n+1<<" j= "<<n<<"\n";
		for(j=1;j<=n;j++){
			for(i=0;i<=n;i++){
				f1<<node[i][j].U_x<<"\t"<<node[i][j].U_y<<"\t"<<"2"<<"\n";
				}
			}
	f1<<"zone "<<"i= "<<n<<" j= "<<n+1<<"\n";
		for(j=0;j<=n;j++){
			for(i=1;i<=n;i++){
				f1<<node[i][j].V_x<<"\t"<<node[i][j].V_y<<"\t"<<"3"<<"\n";
				}
			}
}

void residual_print(){
	s_U = 0.0;
	for(j=1;j<=n;j++){
		for(i=1;i<=n-1;i++){
			s_U = s_U + pow(node[i][j].U - node[i][j].U_old,2);
			}
		}
	s_V =0.0;
	for(j=1;j<=n-1;j++){
		for(i=1;i<=n;i++){
			s_V = s_V + pow(node[i][j].V - node[i][j].V_old,2);
		}
	}
	s_T =0.0;
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){
			s_T = s_T + pow(node[i][j].T - node[i][j].T_old,2);
		}
	}
}

void P_tilda_store(){
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){
			node[i][j].P_tilda_old = node[i][j].P_tilda;
		}
	}
}
void U_plot(){
	ofstream f3;
	f3.open("u_plot_10^6.plt");
	f3<<"zone "<<"i="<<n+1<<" j="<<n<<"\n";
	for(j=1;j<=n;j++){
		for(i=0;i<=n;i++){
			f3<<node[i][j].U_x<<"\t"<<node[i][j].U_y<<"\t"<<node[i][j].U<<"\n";
		}

	}
}
void V_plot(){
	ofstream f4;
	f4.open("V_plot_10^6.plt");
	f4<<"zone "<<"i="<<n<<" j="<<n+1<<"\n";
	for(j=0;j<=n;j++){
		for(i=1;i<=n;i++){
			f4<<node[i][j].V_x<<"\t"<<node[i][j].V_y<<"\t"<<node[i][j].V<<"\n";
		}
	}
}
void T_plot(){
	ofstream f5;
	f5.open("T_plot_10^6.plt");
	f5<<"zone "<<"i="<<n<<" j="<<n<<"\n";
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){
			f5<<node[i][j].P_x<<"\t"<<node[i][j].P_y<<"\t"<<node[i][j].T<<"\n";
		}
	}
}

void CT(){
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){

			node[i][j].uedtedx = node[i][j].U*(node[i][j].T + node[i+1][j].T)/(2.0*dx);
			node[i][j].uwdtwdx = node[i-1][j].U*(node[i][j].T + node[i-1][j].T)/(2.0*dx);
			node[i][j].vndtndy = node[i][j].V*(node[i][j].T + node[i][j+1].T)/(2.0*dy);
			node[i][j].vsdtsdy = node[i][j-1].V*(node[i][j].T + node[i][j-1].T)/(2.0*dy);
			node[i][j].Ct = node[i][j].uedtedx - node[i][j].uwdtwdx + node[i][j].vndtndy - node[i][j].vsdtsdy;
		}
	}
}

void DT(){
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){

			node[i][j].d2tedx2 = (node[i+1][j].T - node[i][j].T)/(dx*dx);
			node[i][j].d2twdx2 = (node[i][j].T - node[i-1][j].T)/(dx*dx);
			node[i][j].d2tndy2 = (node[i][j+1].T - node[i][j].T)/(dy*dy);
			node[i][j].d2tsdy2 = (node[i][j].T - node[i][j-1].T)/(dy*dy);

			node[i][j].Dt = node[i][j].d2tedx2 - node[i][j].d2twdx2 + node[i][j].d2tndy2 - node[i][j].d2tsdy2;
		}
	}
}
void temp(){
	temp_boundary_conditions();
	CT();
	DT();
	t_t = 0.0;
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){
			node[i][j].T_new = node[i][j].T - dt*node[i][j].Ct + dt*node[i][j].Dt;
			t_t = t_t + pow(node[i][j].T_new - node[i][j].T,2);
			node[i][j].T = node[i][j].T_new;
		}
	}
}

void temp_boundary_conditions(){

	//west wall
	for(j=1;j<=n;j++){
		node[0][j].T = 2.0 - node[1][j].T;
	}
	//east wall
	for(j=1;j<=n;j++){
		node[n+1][j].T = -node[n][j].T;
	}
	//north wall
	for(i=1;i<=n;i++){
		node[i][n+1].T = node[i][n].T;
	}
	//south wall
	for(i=1;i<=n;i++){
		node[i][0].T = node[i][1].T;
	}
}
void vector_plot(){
	ofstream f6;
	f6.open("vector_plot_10^6.plt");
	f6<<"zone "<<"i="<<n<<" j="<<n<<"\n";
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++){
			f6<<node[i][j].P_x<<"\t"<<node[i][j].P_y<<"\t"<<((node[i][j].U + node[i-1][j].U)/2.0)<<"\t"<<((node[i][j].V + node[i][j-1].V)/2.0)<<"\n";
		}
	}
}
void Nusselt_plot(){
	//west wall
	ofstream f9;
	f9.open("average number_10^6.txt");
	for(j=1;j<=n;j++){
		node[1][j].Nu = (node[0][j].T-node[2][j].T)/(2*dx);
		node[n][j].Nu = (node[n+1][j].T-node[n-1][j].T)/(2*dx);
		nu_wall_west = nu_wall_west + node[1][j].Nu*dx;
		nu_wall_east = nu_wall_east + node[n][j].Nu*dx;
	}
	cout<<"west wall ="<<nu_wall_west<<" east wall ="<<nu_wall_east;
	f9<<"outerloop"<<outer_loop<<"\t"<<"west wall ="<<nu_wall_west<<" east wall ="<<nu_wall_east;


	ofstream f7,f8;
	f7.open("west_wall_10^4.plt");
	f8.open("east_wall_10^4.plt");

	for(j=1;j<=n;j++){
		f7<<node[1][j].P_y<<"\t"<<node[1][j].Nu<<"\n";
		f8<<node[1][j].P_y<<"\t"<<node[n][j].Nu<<"\n";
	}
}
