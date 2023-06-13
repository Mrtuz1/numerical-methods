#include <stdio.h>
#include <math.h>
#define N 50

void bisection();
void regulaFalsi();
void NewtonRapshon();
void inverseMatrix();
void GaussElemination();
void Gauss_seidel();
void SayisalTurev();
void trapez();
void simpson();
void GregoryNewtonEnterpolasyon();

//-----------------------------------------//

int fonkDer();
void fonksiyonAl(float fonksiyon[N],int derece);
void matrixAl(float matrix[N][N],int n);
double fonkHesap(int derece,float fonksiyon[N],double value);
void sinirAl(double *alts,double *usts);
double turevHesap(int derece,float fonksiyon[N],double value);
float detHesap(float matrix[N][N],int n);
void denkMatrixAl(float denkMatrix[N][N],int denklems,int degiskens);
int factorial(int value);
float C(int n,int r);

int main() {
	int option;
	do {
		printf("\nyontem seciniz\n\n0. cikis\n1. Bisection yontemi\n2. Regula-Falsi yontemi\n3. Newton-Rapshon yontem\n4. NxN lik bir matrisin tersi\n5. Gauss eliminasyon yontemi\n6. Gauss-Seidel yontemi\n7. Sayisal Turev\n8. Simpson yontemi\n9. Trapez yontem\n10. Degisken donusumsuz Gregory-Newton enterpolasyonu\n");
		scanf("%d",&option);
		if(option==1) {
			bisection();
		}
		else if(option==2) {
		    regulaFalsi();
		}
		else if(option==3) {
			NewtonRapshon();
		}
		else if(option==4) {
			inverseMatrix();
		}
		else if(option==5) {
			GaussElemination();
		}
		else if(option==6){
			Gauss_seidel();
		}
		else if(option==7){
			SayisalTurev();
		}
		else if(option==8) {
			simpson();
		}
		else if(option==9) {
			trapez();
		}
		else if(option==10){
			GregoryNewtonEnterpolasyon();
		}

	} while(option!=0);


	return 0;
}


void bisection(){
	double alts,usts,hata,orta;
	int derece=fonkDer(),counter=0;
	float fonksiyon[N];
	fonksiyonAl(fonksiyon,derece);
	sinirAl(&alts,&usts);
	
	if(fonkHesap(derece,fonksiyon,alts)==0 && fonkHesap(derece,fonksiyon,usts)==0) {
		printf("%lf ve %lf fonksiyonun kokleridir!\n",alts,usts);
	} else if(fonkHesap(derece,fonksiyon,alts)==0) {
		printf("%lf fonksiyonun kokudur!\n",alts);
	} else if(fonkHesap(derece,fonksiyon,usts)==0) {
		printf("%lf fonksiyonun kokudur!\n",usts);
	} else if((fonkHesap(derece,fonksiyon,alts)*fonkHesap(derece,fonksiyon,usts))<0) {
		printf("fonksiyonun sinirleri arasinda koku vardir\n");
	} else {
		printf("fonksiyonun sinirleri arasinda koku yoktur ya da sinirleri yanlis girdiniz!\n");
	}
	printf("hatayi giriniz\n");
	if((fonkHesap(derece,fonksiyon,alts)*fonkHesap(derece,fonksiyon,usts))<0){
		
	scanf("%lf",&hata);
		do{
			orta = (alts+usts)/2;
			if(fonkHesap(derece,fonksiyon,orta)*fonkHesap(derece,fonksiyon,alts)<=0){
				usts=orta;
			}
			if(fonkHesap(derece,fonksiyon,orta)*fonkHesap(derece,fonksiyon,usts)<=0){
				alts=orta;
			}
			counter++;
			printf("\n%d. iterasyon\nalts=%f\nusts=%f\norta=%f\n",counter,alts,usts,orta);
		}while(fabs(fonkHesap(derece,fonksiyon,orta))>hata);
		printf("\n%lf fonksiyonun bir kokudur\n",orta);
	}
}


void regulaFalsi(){
	double alts,usts,hata,orta;
	int derece=fonkDer(),counter=0;
	float fonksiyon[N];
	fonksiyonAl(fonksiyon,derece);
	sinirAl(&alts,&usts);
	
	if(fonkHesap(derece,fonksiyon,alts)==0 && fonkHesap(derece,fonksiyon,usts)==0) {
		printf("%lf ve %lf fonksiyonun kokleridir!\n",alts,usts);
	} else if(fonkHesap(derece,fonksiyon,alts)==0) {
		printf("%lf fonksiyonun kokudur!\n",alts);
	} else if(fonkHesap(derece,fonksiyon,usts)==0) {
		printf("%lf fonksiyonun kokudur!\n",usts);
	} else if((fonkHesap(derece,fonksiyon,alts)*fonkHesap(derece,fonksiyon,usts))<0) {
		printf("fonksiyonun sinirleri arasinda koku vardir\n");
	} else {
		printf("fonksiyonun sinirleri arasinda koku yoktur ya da sinirleri yanlis girdiniz!\n");
	}
	if((fonkHesap(derece,fonksiyon,alts)*fonkHesap(derece,fonksiyon,usts))<0){
	printf("hatayi giriniz\n");
	
	scanf("%lf",&hata);
	do{
		orta = ((fonkHesap(derece,fonksiyon,alts)*usts)-(fonkHesap(derece,fonksiyon,usts)*alts))/(fonkHesap(derece,fonksiyon,alts)-fonkHesap(derece,fonksiyon,usts));
		if(fonkHesap(derece,fonksiyon,orta)*fonkHesap(derece,fonksiyon,alts)<=0){
			usts=orta;
		}
		if(fonkHesap(derece,fonksiyon,orta)*fonkHesap(derece,fonksiyon,usts)<=0){
			alts=orta;		
			}
			counter++;
			printf("\n%d. iterasyon\nalts=%f\nusts=%f\norta=%f\n",counter,alts,usts,orta);
		}while(fabs(fonkHesap(derece,fonksiyon,orta))>hata);
		printf("%lf fonksiyonun bir kokudur\n",orta);
	}
}


void NewtonRapshon(){
	double alts,usts,hata,tempx,a,b;;
    int derece=fonkDer(),counter=0;
	float fonksiyon[N];
	fonksiyonAl(fonksiyon,derece);
	sinirAl(&alts,&usts);
	
	printf("baslangic degerini girmek isterseniz 1 e basiniz\n");
	double x;
	scanf("%lf",&x);
	if(x==1){
		printf("baslangýc noktasini giriniz\n");
	    scanf("%lf",&x);
	}
	else{
		x=alts;
	}
	printf("hatayi giriniz\n");
	scanf("%lf",&hata);
	do{
		counter++;
		tempx=x;
		a=fonkHesap(derece,fonksiyon,tempx);
		b=turevHesap(derece,fonksiyon,tempx);
		x=tempx-(a/b);
		
		printf("\n%d. iterasyon\neski x degeri=%f\nyeni x degeri=%f\nf(%f)=%f\nf'(%f)=%f\n",counter,tempx,x,tempx,fonkHesap(derece,fonksiyon,tempx),tempx,turevHesap(derece,fonksiyon,tempx));
		
	}while(fabs(x-tempx)>hata && alts<x && usts>x);
	if(fabs(x-tempx)<hata){
		printf("%lf fonksiyonun bir kokudur\n",x);
	}
	else{
		printf("iraksadi!\n");
	}
}


void inverseMatrix(){
	int n,x,y;
	int satirI,sutunI;
	float determinant;
	float inverseM[N][N];
	float matrix[N][N];
	float kofaktor[N][N];
	float ekmatrix[N][N];
	int i,j,k,l;
	printf("matsisin boyutunu giriniz\n");
	scanf("%d",&n);

	matrixAl(matrix,n);
	determinant=detHesap(matrix,n);
	printf("determinant=%f\n",determinant);
	if(determinant!=0){
		
		for(i=0;i<n;i++){
			
			for(j=0;j<n;j++){
				
				satirI=0;
				for(k=0;k<n;k++){
					sutunI=0;
					for(l=0;l<n;l++){
						
						if(k!=i && l!=j){
							kofaktor[satirI][sutunI]=matrix[k][l];
							sutunI++;
						}
					}
					if(sutunI==n-1){
						satirI++;
					}
				}
				if((i+j)%2==0){
					ekmatrix[i][j]=detHesap(kofaktor,n-1)/determinant;
				}
				else{
					ekmatrix[i][j]=-detHesap(kofaktor,n-1)/determinant;
				}
				
			}
		}
		
		printf("matrisin tersi=\n");
		for(i=0;i<n;i++){
			printf("\n");
			for(j=0;j<n;j++){
				inverseM[i][j]=ekmatrix[j][i];
				printf(" %lf ",inverseM[i][j]);
			}
		}
		
	}
	else{
		printf("matrisin tersi yoktur ! determinant=%lf\n",determinant);
	}
}


void GaussElemination(){
	float x[N];
	int i,j,k,counter=0,X,Y,counter1,coor[N],n;
	float denkMatrix[N][N],temp,maks;
	int denklems,degiskens;
	float kats,a;
	
	
	printf("denklem ve degisken sayisi olan n sayisini giriniz\n");
	scanf("%d",&n);
	
    
	denkMatrixAl(denkMatrix,n,n);
	printf("kosegenlerin carpiminin en buyuk oldugu sekilde matris duzenlenir\n");
	counter=0;	
	do{
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			counter1=1;
			for(k=0;k<counter;k++){
				if(i==coor[k] || j==coor[k]){
					counter1=0;
				}
				
			}
			if(counter1==1 && maks<fabs(denkMatrix[i][j])){
				maks=fabs(denkMatrix[i][j]);
				X=i;
				Y=j;
				coor[counter]=j;
				
			}
		}
	}
	
	maks=0;
	for(i=0;i<=n;i++){
		temp=denkMatrix[Y][i];
		denkMatrix[Y][i]=denkMatrix[X][i];
		denkMatrix[X][i]=temp;	
	}
	counter++;
	
	for(i=0;i<n;i++){
		for(j=0;j<=n;j++){
			printf("%f ",denkMatrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	
	
}while(counter<n);
    
    for(i=0;i<n-1;i++){
    	kats=denkMatrix[i][i];
    	
	for(j=i+1;j<n;j++){
		a=denkMatrix[j][i]/kats;
		for(k=i;k<n+1;k++){
			denkMatrix[j][k]-=denkMatrix[i][k]*a;
		}
	}
}
    if(denkMatrix[n-1][n-1]==0){
    	for(i=0;i<n;i++){
    		x[i]=0;
		}
	}
	else{
		counter=0;
		for(i=n-1;i>=0;i--){
			counter++;
			for(j=n-1;j>n-counter;j--){
				denkMatrix[i][n]-=denkMatrix[i][j]*x[j];
			}
			x[i]=denkMatrix[i][n]/denkMatrix[i][i];
		}
		for(i=0;i<n;i++){
			printf("%d.degisken = %f\n",i+1,x[i]);
		}
	}
}


void Gauss_seidel(){
	int i,j,k,x,y,x1,y1;
	float maks=0,sum;
	float xdeger[2][N],hata;
	float denkMatrix[N][N];
	int n,counter,counter1,iterasyon=0;
	float temp;
	int coor[N];
	
	printf("denklem sayisi ve degisken sayisi olan n sayisini girinz!\n");
	scanf("%d",&n);
    
	denkMatrixAl(denkMatrix,n,n);
	printf("kosegenlerinin carpimi en buyuk olacak sekilde satirlar siralanir\n");
	counter=0;	
	do{
		
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			counter1=1;
			for(k=0;k<counter;k++){
				if(i==coor[k] || j==coor[k]){
					counter1=0;
				}
				
			}
			if(counter1==1 && maks<fabs(denkMatrix[i][j])){
				maks=fabs(denkMatrix[i][j]);
				x=i;
				y=j;
				coor[counter]=j;
				
			}
		}
	}
	maks=0;
	for(i=0;i<=n;i++){
		temp=denkMatrix[y][i];
		denkMatrix[y][i]=denkMatrix[x][i];
		denkMatrix[x][i]=temp;	
	}
	counter++;
}while(counter<n);


for(i=0;i<n;i++){
		for(j=0;j<=n;j++){
			printf("%f ",denkMatrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");

		printf("hata degerini giriniz\n");
		scanf("%f",&hata);
		printf("degiskenlerin baslangic degerlerini giriniz\n");
		for(i=0;i<n;i++){
			printf("%d. degiskenin degerini giriniz\n",i+1);
			scanf("%f",&xdeger[0][i]);
		}
		
		do{
			for(i=0;i<n;i++){
				sum=denkMatrix[i][n];
				for(j=0;j<n;j++){
					if(j!=i){
					sum-=denkMatrix[i][j]*xdeger[0][j];
					}
					
				}
				temp=xdeger[0][i];
				xdeger[0][i]=sum/denkMatrix[i][i];
				xdeger[1][i]=temp;
			}
			
			counter=0;
			for(i=0;i<n;i++){
				if(fabs(xdeger[0][i]-xdeger[1][i])>hata){
					counter=1;
				}
			}
			iterasyon++;
			printf("%d.iterasyon\n",iterasyon);
			for(i=0;i<n;i++){
		    printf("%d. degiskenin degeri = %f\n",i+1,xdeger[0][i]);
	        }
			
		
			
		}while(counter==1);
		
	
}

void simpson(){
	float n,h,i,S=0,fonksiyon[N],x,counter,a=0,b=0,exs,x1,x2,temp;
	int opt;
	double alts,usts;
	int derece=fonkDer();
	fonksiyonAl(fonksiyon,derece);
	sinirAl(&alts,&usts);
	
	
	printf("simpson 1/3=\n");
	
	printf("n degeri giriniz(n cift olmali!)\n");
	scanf("%f",&n);
	h=(usts-alts)/n;
	printf("h=%f\n",h);
	i=alts+h;
	while(i<usts){
		S+=4*fonkHesap(derece,fonksiyon,i);
		i+=2*h;
	} 
	i=alts+2*h;
	while(i<usts){
		S+=2*fonkHesap(derece,fonksiyon,i);
		i+=2*h;
	}
	S+=fonkHesap(derece,fonksiyon,alts)+fonkHesap(derece,fonksiyon,usts);
	S=(S*h/3);
	printf("simpson 1/3 integral= %f\n\n",S);
	S=0;
	printf("simpson 3/8=\n");
		printf("n degerini giriniz\n");
		scanf("%f",&n);
		h=(usts-alts)/3;
		
		x=(fabs(usts)+fabs(alts))/n;
		counter=0;
		if(n==1){
			S=(usts-alts)*(fonkHesap(derece,fonksiyon,alts)+fonkHesap(derece,fonksiyon,usts)+3*(fonkHesap(derece,fonksiyon,alts+h)+fonkHesap(derece,fonksiyon,alts+2*h)))/8;
		}
		else{
			do{
				if(counter==0){
					a=alts;
					b=a+x;
					x1=a+(b-a)/3;
					x2=a+(b-a)*2/3;
					S+=(b-a)*(fonkHesap(derece,fonksiyon,a)+fonkHesap(derece,fonksiyon,b)+3*(fonkHesap(derece,fonksiyon,x1)+fonkHesap(derece,fonksiyon,x2)))/8;
					exs=b;
						printf("alt sinir=%f  ust sinir=%f  integralin degeri=%f\n",a,b,((b-a)*(fonkHesap(derece,fonksiyon,a)+fonkHesap(derece,fonksiyon,b)+3*(fonkHesap(derece,fonksiyon,x1)+fonkHesap(derece,fonksiyon,x2)))/8));
				}
				else if(counter==n-1){
					b=usts;
					a=b-x;
					x1=a+(b-a)/3;
					x2=a+(b-a)*2/3;
					S+=(b-a)*(fonkHesap(derece,fonksiyon,a)+fonkHesap(derece,fonksiyon,b)+3*(fonkHesap(derece,fonksiyon,x1)+fonkHesap(derece,fonksiyon,x2)))/8;
						printf("alt sinir=%f  ust sinir=%f  integralin degeri=%f\n",a,b,((b-a)*(fonkHesap(derece,fonksiyon,a)+fonkHesap(derece,fonksiyon,b)+3*(fonkHesap(derece,fonksiyon,x1)+fonkHesap(derece,fonksiyon,x2)))/8));
				}
				else{
					a=exs;
					b=a+x;
					x1=a+(b-a)/3;
					x2=a+(b-a)*2/3;
					S+=(b-a)*(fonkHesap(derece,fonksiyon,a)+fonkHesap(derece,fonksiyon,b)+3*(fonkHesap(derece,fonksiyon,x1)+fonkHesap(derece,fonksiyon,x2)))/8;
					exs=b;
						printf("alt sinir=%f  ust sinir=%f  integralin degeri=%f\n",a,b,((b-a)*(fonkHesap(derece,fonksiyon,a)+fonkHesap(derece,fonksiyon,b)+3*(fonkHesap(derece,fonksiyon,x1)+fonkHesap(derece,fonksiyon,x2)))/8));
				}
				counter++;
			}while(counter<n);
		}
		printf(" simpson 3/8 integral=%f\n",(S));
}



void SayisalTurev(){
	int derece;
	float x,h,turev;
	float fonksiyon[N];
	printf("fonksiyonun derecesini giriniz\n");
	scanf("%d",&derece);
	fonksiyonAl(fonksiyon,derece);
	printf("turevini istediginiz x degerini giriniz\n");
	scanf("%f",&x);
	printf("h degeri giriniz(epsilon)\n");
	scanf("%f",&h);
	
	

	turev=(fonkHesap(derece,fonksiyon,x)-fonkHesap(derece,fonksiyon,x-h))/h;
	printf("geri fark= %f\n",turev);
	
	turev=(fonkHesap(derece,fonksiyon,x+h)-fonkHesap(derece,fonksiyon,x))/h;
	printf("ileri fark= %f\n",turev);
	
	turev=(fonkHesap(derece,fonksiyon,x+h)-fonkHesap(derece,fonksiyon,x-h))/(2*h);
	printf("merkezi fark= %f\n",turev);
	
}



void GregoryNewtonEnterpolasyon(){
	
	int n,i,j;
    float X[N],Y[N],sum,h,a,x;
    float fark[N][N];
    
    printf("\n(kac tane x var?)n degerini giriniz:\n");
    scanf("%d",&n);
    printf("x ve y degerlerini sirasiyla giriniz!\n");
    for(i=0;i<n;i++){
    	printf("%d. x ve y degerlerini gir\n",i+1);
        scanf("%f",&X[i]);
        scanf("%f",&Y[i]);
    }
    
    printf("sonucunu istediginiz x degerini giriniz\n");
    scanf("%f",&x);
    
    //fark tablosu
    for(i=0;i<n-1;i++){
    	fark[i][0]=Y[i+1]-Y[i];
	}
	for(i=1;i<4;i++){
		for(j=0;j<n-i-1;j++){
			fark[j][i]=fark[j+1][i-1]-fark[j][i-1];
			
		}
	}
	for(i=0;i<n;i++){
		printf("%f-->%f\n",X[i],Y[i]);
	}
	printf("\nfark tablosu=\n");
	for(i=0;i<n-1;i++){
		for(j=0;j<n-i-1;j++){
			printf("%f ",fark[i][j]);
		}
		printf("\n");
	}
	
	
	
    h=X[1]-X[0];
    printf("h=&f",h);
    
    sum=Y[0];
    a=1;
    for(i=0;i<4;i++){
    	a=(a*(x-X[i]))/(h*(i+1));
    	sum+=a*fark[0][i];
    	
	}
	printf("sonuc= %f",sum);
}


void trapez(){
	float n,h,i,S,fonksiyon[N];
	double alts,usts;
	int derece=fonkDer();
	fonksiyonAl(fonksiyon,derece);
	sinirAl(&alts,&usts);
	printf("n degeri giriniz\n");
	scanf("%f",&n);
	h=(usts-alts)/n;
	i=alts+h;
	while(i<usts){
			S+=fonkHesap(derece,fonksiyon,i);
			i+=h;
	}
	S+=(fonkHesap(derece,fonksiyon,alts)+fonkHesap(derece,fonksiyon,usts))/2;
	S=(S*h);
	printf("integral=%f\n",S);
	
}


int fonkDer(){
	int derece;
	printf("\nFonksiyonunuzun derecesi kac?\n");
	scanf("%d", &derece);
	return derece;
}

void fonksiyonAl(float fonksiyon[N],int derece){
	printf("fonksiyonu giriniz\n");
	int i;	
	for(i=0;i<=derece;i++){
		printf("%d. elemanin katsayisini giriniz?\n", i);
		scanf("%f",&fonksiyon[i]);
	}
	printf("fonksiyon=");
	for(i=0;i<=derece;i++){
		printf("(%f*x^%d)",fonksiyon[i],i);
		if(i!=derece){
			printf("+");
		}
		
	}
	printf("\n");
	
}

double fonkHesap(int derece,float fonksiyon[N],double value) {     //x=FonkHesap(derece,fonksiyon,3.8);
	int i;
	double sonuc=0.0;
	for(i=0;i<=derece;i++){
		sonuc+=fonksiyon[i]*pow(value,i);
	}
	return sonuc;
}

void sinirAl(double *alts,double *usts){
 	double tempA,tempU;
 	do{
 	printf("alt siniri giriniz\n");
 	scanf("%lf",&tempA);
 	printf("ust siniri giriniz\n");
 	scanf("%lf",&tempU);
 	if(tempA>tempU){
 	printf("alt sinir ust sinirdan buyuk olamaz!\n");
	}
	}while(tempA>tempU);
 	*alts=tempA;
 	*usts=tempU;
 	
}

double turevHesap(int derece,float fonksiyon[N],double value){
	int i;
	double sonuc=0.0;
	for(i=1;i<=derece;i++){
		sonuc+=fonksiyon[i]*i*pow(value,i-1);
	}
	
	return sonuc;
}

void matrixAl(float matrix[N][N],int n){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			printf("%d.satir %d.sutundaki elemani giriniz\n",i+1,j+1);
			scanf("%f",&matrix[i][j]);
		}
	}
	printf("matrix=\n");
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			printf("%f ",matrix[i][j]);
		}
		printf("\n");
	}
	
}

float detHesap(float matrix[N][N],int n){
	float x;
	float minor[N][N];
	float det;
	int i,j,k;
	int satirI,sutunI;
	if(n==1){
		det=matrix[0][0];
	}
	else if(n==2){
		det=(matrix[1][1]*matrix[0][0])-(matrix[1][0]*matrix[0][1]);
	}
	else if(n>2){
		det=0;
		for(i=0;i<n;i++){// i ilk satýrdaki index olacak
		satirI=0;
		//matrix[0][i] nin minorunu alacak!
		    for(j=1;j<n;j++){// satýrlarý tarayacak
		    sutunI=0;
		    	for(k=0;k<n;k++){//sütunlarý tarayacak
		    		if(k!=i){
		    			minor[satirI][sutunI]=matrix[j][k];
		    			sutunI++;
					}
				}
				satirI++;
			}
			
			x=detHesap(minor,n-1);

			if(i%2==0){
				det+=matrix[0][i]*x;
			}
			else{
				det-=matrix[0][i]*x;
			}	

			
		}
	}
	
return det;
}



void denkMatrixAl(float denkMatrix[N][N],int denklems,int degiskens){
	
	int i,j;
	
	for(i=0;i<denklems;i++){
		for(j=0;j<degiskens;j++){
			printf("%d.denklemin %d.degiskeninin katsayisini giriniz!\n",i+1,j+1);
			scanf("%f",&denkMatrix[i][j]);
		}
	}
	
	for(i=0;i<denklems;i++){
		printf("%d.denklemin sonucunu yaziniz!\n",i+1);
		scanf("%f",&denkMatrix[i][degiskens]);
	}
}

int factorial(int value){
	int i,sum=1;
	for (i=1;i<=value;i++){
		
		sum=sum*i;
		
	}
	
	return sum;
	
}

float C(int n,int r){
	
	float sum=1;
	int i;
	
	for(i=n;i>r;i--){
		sum=sum*i;
	}
	sum = sum /factorial(r);
	
	
	return sum;
	
}
