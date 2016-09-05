import java.util.*;

public class CorrelationLeaveOneOut {

 public static void main(String args[]){
		float x1[] = {86  ,  97  ,  99  , 100 ,  101 ,  103 ,  106 ,  110 ,  112 ,  113};
		float x2[] = {500,   163,   197,   195,   217,   227,   237,   198,   242,   242};
		System.out.println("Pearson corr = "+calcCorrelationCoeff(x1, x2));
		System.out.println("Spearman corr = "+calcSpearmanCorrelationCoeff(x1, x2));
		System.out.println("Robust Pearson corr = "+ calcCorrelationCoeffLeaveOneOut(x1, x2, 2f, 0));
		System.out.println("Robust Spearman corr = "+ calcCorrelationCoeffLeaveOneOut(x1, x2, 2f, 1));
 }
 
 /**
  * 
  * @param m1
  * @param m2
  * @param ZValueThreshold
  * @param corrType  0 - for Pearson correlation, 1 - for Spearman correlation
  * @return
  */
 public static float calcCorrelationCoeffLeaveOneOut(float m1[], float m2[], float ZValueThreshold, int corrType){
	  float res = 0;
	  
	  float correlationCoeffs[] = new float[m1.length];

	  float m1s[] = new float[m1.length-1];
	  float m2s[] = new float[m2.length-1];
	  
	  for(int i=0;i<m1.length;i++){
		  int k=0;
		  for(int j=0;j<m1.length;j++){
			  if(j!=i){
				  m1s[k] = m1[j];
				  m2s[k] = m2[j];
				  k++;
			  }
		  }
		  float corrCoeff = 0f;
		  switch(corrType){
		  	case 0: corrCoeff = calcCorrelationCoeff(m1s,m2s); break;
		  	case 1: corrCoeff = calcSpearmanCorrelationCoeff(m1s,m2s); break;
		  }
		  correlationCoeffs[i] = corrCoeff;
	  }
	  
	  Vector m1new = new Vector();
	  Vector m2new = new Vector();
	  
	  float corrCoeff1[] = new float[m1.length-1];
	  for(int i=0;i<m1.length;i++){
		  int k=0;
		  for(int j=0;j<m1.length;j++){
			  if(j!=i){
				  corrCoeff1[k] = correlationCoeffs[j];
				  k++;
			  }
		  }
		  float meanValue = calcMean(corrCoeff1);
		  float stdDev = calcStandardDeviation(corrCoeff1);
		  float zvalue = (correlationCoeffs[i]-meanValue)/stdDev;
		  if(Math.abs(zvalue)<ZValueThreshold){
			 m1new.add(new Float(m1[i])); 
			 m2new.add(new Float(m2[i]));
		  }
	  }
	  
	  m1s = new float[m1new.size()];
	  m2s = new float[m2new.size()];
	  for(int i=0;i<m1new.size();i++){
		  m1s[i] = ((Float)m1new.get(i)).floatValue();
		  m2s[i] = ((Float)m2new.get(i)).floatValue();
	  }

	  switch(corrType){
	  	case 0: res = calcCorrelationCoeff(m1s,m2s); break;
	  	case 1: res = calcSpearmanCorrelationCoeff(m1s,m2s); break;
	  }
	  
	  /*## Function to compute robust correlation coeffients between mRNA-microRNA expression profiles
	  ## Input parameters:
	  ## a,b = vectors values 
	  ## lth = P-value threshold for robust correlation
	  lout <- function(x,y,lth)
	  {
	  	x <- as.numeric(x)
	  	y <- as.numeric(y)
	  	lth <- lth 
	  	pear.out <- c()
	  	spe.out <- c()
	  	for (k in 1:length(x))
	  	{
	  		tempx <- c()
	  		tempy <- c()
	  		for (z in 1:length(x))
	  		{
	  			if (z != k)
	  			{
	  				tempx <- c(tempx, x[z])
	  				tempy <- c(tempy, y[z])
	  			}
	  		}
	  		testpear <- cor.test(tempx,tempy, method = "pearson")
	  		pear.out <- c(pear.out, as.numeric(testpear$estimate))
	  		testspe <- cor.test(tempx,tempy, method = "spearman")
	  		spe.out <- c(spe.out, as.numeric(testspe$estimate))
	  	}
	  	newx <- c()
	  	newy <- c()
	  	for (k in 1:length(x))
	  	{
	  		ptemp <- pnorm(pear.out[k], mean(pear.out), sd(pear.out), lower.tail = TRUE, log.p = FALSE)
	  		pv <- 2*(1 - ptemp)
	  		if (pv >= lth)
	  		{
	  			newx <- c(newx, x[k])
	  			newy <- c(newy, y[k])
	  		}
	  	}
	  	pear.tot <<- cor.test(x,y, method = "pearson")
	  	spe.tot <<- cor.test(x,y, method = "spearman")
	  	pear.rob <<- cor.test(newx,newy, method = "pearson")
	  	spe.rob <<- cor.test(newx,newy, method = "spearman")
	  	return("ok")
	  }*/
	  
	  
	  return res;
 }
 
	
	
  public static float calcCorrelationCoeff(float m1[], float m2[]){
  float res = 0;
  int N = m1.length;
  float xy = 0f, x2 = 0f, y2 = 0f, x = 0f, y = 0f;
  for(int i=0;i<N;i++){
    xy+=m1[i]*m2[i];
    x+=m1[i];
    y+=m2[i];
    x2+=m1[i]*m1[i];
    y2+=m2[i]*m2[i];
  }
  double disp = Math.sqrt((x2-x*x/N)*(y2-y*y/N));
  if(Math.abs(disp)<1e-20)
	  res = 0;
  else
      res = (float)((xy-x*y/N)/disp);
  return res;
  }
  
  public static double calcCorrelationPValue(float correlationValue, int numberOfPoints){
	  double pvalue = 0;
	  int degreeOfFreedom = numberOfPoints-2;
	  if(Math.abs(correlationValue)<1)if(degreeOfFreedom>2){
		  double tt = correlationValue/Math.sqrt((1-correlationValue*correlationValue)/(degreeOfFreedom));
		  //System.out.println("corr="+correlationValue+" tt="+tt);
		  if(Math.abs(correlationValue)>=1) pvalue=0;		  
		  if(degreeOfFreedom==2) pvalue=1;
		  pvalue = (1-tcdf(tt,degreeOfFreedom))*2;
	  }
	  return pvalue;
  }
  


  public static float calcSpearmanCorrelationCoeff(float m1[], float m2[]){
  float res = 0;
  int N = m1.length;
  int ind1[] = SortMass(m1);
  int ind2[] = SortMass(m2);
  
  int rank1[] = ind2rank(ind1);  
  int rank2[] = ind2rank(ind2);
  
  //for(int i=0;i<ind2.length;i++)
  //   System.out.println(rank2[i]);
  
  float rank1f[] = new float[ind1.length];
  float rank2f[] = new float[ind1.length];
  for(int i=0;i<ind1.length;i++) rank1f[i] = (float)rank1[i];
  for(int i=0;i<ind2.length;i++) rank2f[i] = (float)rank2[i];
  /*long d2 = 0;
  for(int i=0;i<ind1.length;i++){
    d2+=(ind1[i]-ind2[i])*(ind1[i]-ind2[i]);
  }
  res = 1f-6f*(float)d2/((float)N*((float)N*(float)N-1f));*/
  res = calcCorrelationCoeff(rank1f,rank2f);
  return res;
  }
  
  public static int[] ind2rank(int ind[]){
	  int rank[] = new int[ind.length];
	  for(int i=0;i<ind.length;i++)
		  rank[ind[i]] = i;
	  return rank;
  }
  
  public static int[] SortMass(Float mass[]){
	  float m[] = new float[mass.length];
	  for(int i=0;i<mass.length;i++)
		  m[i] = mass[i];
	  return SortMass(m);
  }
  

  public static int[] SortMass(float cais[]){
  int res[]=new int[cais.length];
  for (int i = 0; i < res.length; i++) res[i]=i;

  int i,j,k,inc,n=cais.length;
  float v;

  inc=1;
  do {
  	inc *= 3;
  	inc++;
  } while (inc <= n);

  do {
  	inc /= 3;
  	for (i=inc+1;i<=n;i++) {
  		v=cais[res[i-1]];
  		j=i;
                  k=res[i-1];
  		while (cais[res[j-inc-1]]>v) {
  			//cais[j]=cais[j-inc];
                          res[j-1]=res[j-inc-1];
  			j -= inc;
  			if (j <= inc) break;
  		}
  		//cais[j]=v;
                  res[j-1]=k;
  	}
  } while (inc > 0);

  return res;
}
  

  public static float calcMedian(float f[]){
    float r = 0;
    int ind[] = SortMass(f);
    if(f.length!=0){
    if(f.length==2*(int)(0.5f*f.length)){
      int mid1 = (int)(0.5f*f.length)-1;
      int mid2 = (int)(0.5f*f.length);
      r = 0.5f*(f[ind[mid1]]+f[ind[mid2]]);
    }else{
      int mid = (int)(0.5f*f.length);
      r = f[ind[mid]];
    }}
    return r;
  }

  public static float calcStandardDeviation(float f[]){
    float r = 0;
    float x = 0;
    float x2 = 0;
    for(int i=0;i<f.length;i++){
      x+=f[i];
      x2+=f[i]*f[i];
    }
    x/=f.length;
    r = (float)Math.sqrt((x2/f.length-x*x)*(float)f.length/((float)f.length-1));
    return r;
  }
  
  public static float calcDispersion(float f[]){
	    float r = 0;
	    float x = 0;
	    float x2 = 0;
	    for(int i=0;i<f.length;i++){
	      x+=f[i];
	      x2+=f[i]*f[i];
	    }
	    x/=f.length;
	    r = (float)(x2/f.length-x*x);
	    return r;
	  }
  

  public static float calcMean(float f[]){
    float x = 0;
    for(int i=0;i<f.length;i++){
      x+=f[i];
    }
    return x/f.length;
  }

  /**
   * logarithm of gamma function
   * @param xx
   * @return
   */
  public static double gammln(double xx)
  {
  	int j;
  	double x,y,tmp,ser;
  	double cof[]= {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};

  	y=x=xx;
  	tmp=x+5.5;
  	tmp -= (x+0.5)*Math.log(tmp);
  	ser=1.000000000190015;
  	for (j=0;j<6;j++) ser += cof[j]/++y;
  	return -tmp+Math.log(2.5066282746310005*ser/x);
  }
  /**
   * continued fraction used by betai
   */
  public static double betacf(double a, double b, double x)
  {
  	int MAXIT=100;
  	double EPS=1e-10;
  	double FPMIN=Double.MIN_VALUE/EPS;
  	int m,m2;
  	double aa,c,d,del,h,qab,qam,qap;

  	qab=a+b;
  	qap=a+1.0;
  	qam=a-1.0;
  	c=1.0;
  	d=1.0-qab*x/qap;
  	if (Math.abs(d) < FPMIN) d=FPMIN;
  	d=1.0/d;
  	h=d;
  	for (m=1;m<=MAXIT;m++) {
  		m2=2*m;
  		aa=m*(b-m)*x/((qam+m2)*(a+m2));
  		d=1.0+aa*d;
  		if (Math.abs(d) < FPMIN) d=FPMIN;
  		c=1.0+aa/c;
  		if (Math.abs(c) < FPMIN) c=FPMIN;
  		d=1.0/d;
  		h *= d*c;
  		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
  		d=1.0+aa*d;
  		if (Math.abs(d) < FPMIN) d=FPMIN;
  		c=1.0+aa/c;
  		if (Math.abs(c) < FPMIN) c=FPMIN;
  		d=1.0/d;
  		del=d*c;
  		h *= del;
  		if (Math.abs(del-1.0) <= EPS) break;
  	}
  	if (m > MAXIT) System.out.println("ERROR: a or b too big, or MAXIT too small in betacf");
  	return h;
  }
  /**
   * incomplete beta function
   * @param a
   * @param b
   * @param x
   * @return
   */
  public static double betai(double a, double b, double x)
  {
	double bt;

  	if (x < 0.0 || x > 1.0) System.out.println("ERROR: Bad x in routine betai");
  	if (x == 0.0 || x == 1.0) bt=0.0;
  	else
  		bt=Math.exp(gammln(a+b)-gammln(a)-gammln(b)+a*Math.log(x)+b*Math.log(1.0-x));
  	if (x < (a+1.0)/(a+b+2.0))
  		return bt*betacf(a,b,x)/a;
  	else
  		return 1.0-bt*betacf(b,a,1.0-x)/b;
  }
  /**
   * Student's t-test for difference of means
   * @param t
   * @param df
   * @return
   */
  public static double ttest(double t, int df)
  {
  	return betai(0.5*df,0.5,df/(df+t*t));
  }
  /**
   * Matlab's tcdf
   * @param t
   * @param df
   * @return
   */
  public static double tcdf(double t, int df)
  {
  	return 1-ttest(t,df)/2;
  }


}