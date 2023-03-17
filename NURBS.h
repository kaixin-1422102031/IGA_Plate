int Spatch::FindSpan(int n,int p,double u,vector<double> &U)
{

	if(u==U[n+1])return(n);
	int low=p;
	int high=n+1;
	int mid=(high+low)/2;
	while(u<U[mid]||u>=U[mid+1])
	{
		if(u<U[mid]) high=mid;
		else low=mid;
		mid=(high+low)/2;

	}
	return(mid);

}

//求基函数的值 （函数求值）
double Spatch::OneBasisFun_withreturn(int p,int m,vector<double> &U,int i,double u)
{
	double Nip;
	vector<double> N(p+1);
	if((i==0&&u==U[0])||(i==m-p-1&&u==U[m]))
	{
		Nip=1.0;
		return Nip;
	}
	if(u<U[i]||u>=U[i+p+1])
	{
		Nip=0.0;
		return Nip;
	}
	int jj,j;
	for(jj=0;jj<=p;jj++){
		if(u>=U[i+jj]&&u<U[i+jj+1])
			N[jj]=1.0;
		else 
			N[jj]=0.0;
	}
	double saved;
	int k;
	for(k=1;k<=p;k++)
	{
		if(N[0]==0.0) saved =0.0;
		else saved=((u-U[i])*N[0])/(U[i+k]-U[i]);
		for(j=0;j<p-k+1;j++)
		{
			double Uleft=U[i+j+1];
			double Uright=U[i+j+k+1];
			if(N[j+1]==0.0)
			{
				N[j]=saved;saved=0.0;
			}
			else
			{
				double temp=N[j+1]/(Uright-Uleft);
				N[j]=saved +(Uright-u)*temp;
				saved=(u-Uleft)*temp;
			}
		}

	}
	Nip=N[0];
	return Nip;
	//cout<< Nip <<" ";

}

//求单个基函数的导数
double Spatch::DersOneBasisFun_double_First_Ders(int p,int m,vector<double> &U,int i,double u,int n)
{
	double DersOnebasisfun_value;
	if(p==3){
		double Onebasisfun_i_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i,u);
		double Onebasisfun_iplus1_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i+1,u);
		double Uipi=1.0/(U[i+p]-U[i]);
		double Uip1i1=1.0/(U[i+p+1]-U[i+1]);
		if((U[i+p]-U[i])==0)
			Uipi=0.0;
		if((U[i+p+1]-U[i+1])==0)
			Uip1i1=0.0;
		DersOnebasisfun_value=p*(Onebasisfun_i_pminus1*Uipi-Onebasisfun_iplus1_pminus1*Uip1i1);
		if(u==1.0)
		{
			u=0.99999;
			double Onebasisfun_i_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i,u);
			double Onebasisfun_iplus1_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i+1,u);
			double Uipi=1.0/(U[i+p]-U[i]);
			double Uip1i1=1.0/(U[i+p+1]-U[i+1]);
			if((U[i+p]-U[i])==0)
				Uipi=0.0;
			if((U[i+p+1]-U[i+1])==0)
				Uip1i1=0.0;
			double DersOnebasisfun_value_zero=p*(Onebasisfun_i_pminus1*Uipi-Onebasisfun_iplus1_pminus1*Uip1i1);
			DersOnebasisfun_value=-DersOnebasisfun_value_zero;
			//cout<< "Hehehehhe" <<  " " << i << " " << Onebasisfun_i_pminus1  << " " << Uipi  << " " <<  Onebasisfun_iplus1_pminus1  << " " << Uip1i1 << " " << DersOnebasisfun_value_zero << "\n";
		}

	}

	if(p==2){
		double Onebasisfun_i_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i,u);
		double Onebasisfun_iplus1_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i+1,u);
		double Uipi=1.0/(U[i+p]-U[i]);
		double Uip1i1=1.0/(U[i+p+1]-U[i+1]);
		if((U[i+p]-U[i])==0)
			Uipi=0.0;
		if((U[i+p+1]-U[i+1])==0)
			Uip1i1=0.0;
		DersOnebasisfun_value=p*(Onebasisfun_i_pminus1*Uipi-Onebasisfun_iplus1_pminus1*Uip1i1);
		if(u==1.0)
		{
			u=0.99999;
			double Onebasisfun_i_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i,u);
			double Onebasisfun_iplus1_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i+1,u);
			double Uipi=1.0/(U[i+p]-U[i]);
			double Uip1i1=1.0/(U[i+p+1]-U[i+1]);
			if((U[i+p]-U[i])==0)
				Uipi=0.0;
			if((U[i+p+1]-U[i+1])==0)
				Uip1i1=0.0;
			double DersOnebasisfun_value_zero=p*(Onebasisfun_i_pminus1*Uipi-Onebasisfun_iplus1_pminus1*Uip1i1);
			DersOnebasisfun_value=-DersOnebasisfun_value_zero;
			//cout<< "Hehehehhe" <<  " " << i << " " << Onebasisfun_i_pminus1  << " " << Uipi  << " " <<  Onebasisfun_iplus1_pminus1  << " " << Uip1i1 << " " << DersOnebasisfun_value_zero << "\n";
		}

	}
	if(p==1)
	{
		if(i==0)
			DersOnebasisfun_value=-1;
		else{
			DersOnebasisfun_value=1;
		}
	}

		//double Onebasisfun_i_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i,u);
		//double Onebasisfun_iplus1_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i+1,u);
		//double Uipi=1.0/(U[i+p]-U[i]);
		//double Uip1i1=1.0/(U[i+p+1]-U[i+1]);
		//if((U[i+p]-U[i])==0)
		//	Uipi=0.0;
		//if((U[i+p+1]-U[i+1])==0)
		//	Uip1i1=0.0;
		//DersOnebasisfun_value=p*(Onebasisfun_i_pminus1*Uipi-Onebasisfun_iplus1_pminus1*Uip1i1);
		//if(u==1.0)
		//{
		//	u=0.99999;
		//	double Onebasisfun_i_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i,u);
		//	double Onebasisfun_iplus1_pminus1=OneBasisFun_withreturn(p-1,U.size()-1,U,i+1,u);
		//	double Uipi=1.0/(U[i+p]-U[i]);
		//	double Uip1i1=1.0/(U[i+p+1]-U[i+1]);
		//	if((U[i+p]-U[i])==0)
		//		Uipi=0.0;
		//	if((U[i+p+1]-U[i+1])==0)
		//		Uip1i1=0.0;
		//	double DersOnebasisfun_value_zero=p*(Onebasisfun_i_pminus1*Uipi-Onebasisfun_iplus1_pminus1*Uip1i1);
		//	DersOnebasisfun_value=-DersOnebasisfun_value_zero;
		//	//cout<< "Hehehehhe" <<  " " << i << " " << Onebasisfun_i_pminus1  << " " << Uipi  << " " <<  Onebasisfun_iplus1_pminus1  << " " << Uip1i1 << " " << DersOnebasisfun_value_zero << "\n";
		//}

	return DersOnebasisfun_value;

}

void Spatch::Convert_NURBSfileCP_Into_MatrixControlPoints()
{
	int i,j;
	ControlPoints_Matrix_Type_2D.resize(row, vector <conpoints> (col));
	for(i=0;i<row;i++)
		for(j=0;j<col;j++)
		{
			int index_temp=col*i+j;
			ControlPoints_Matrix_Type_2D[i][j].w=ControlPoints[index_temp].w;
			ControlPoints_Matrix_Type_2D[i][j].x=ControlPoints[index_temp].x;
			ControlPoints_Matrix_Type_2D[i][j].y=ControlPoints[index_temp].y;
			ControlPoints_Matrix_Type_2D[i][j].z=ControlPoints[index_temp].z;
		}
}
