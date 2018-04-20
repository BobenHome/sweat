package beeColony;
import java.lang.Math;

public  class beeColony {
	/* Control Parameters of ABC algorithm（ABC算法的控制参数）*/
	//蜂群大小（采蜜蜂数量+观察蜂数量）
	int NP=20; /* The number of colony size (employed bees+onlooker bees)*/
	//食物源大小（可行解的数量）
	int FoodNumber = NP/2; /*The number of food sources equals the half of the colony size*/
	//无法通过limit次试验来改进的食物源将被采蜜蜂遗弃
	int limit = 100;  /*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
	//觅食的周期数{停止标准}
	int maxCycle = 2500; /*The number of cycles for foraging {a stopping criteria}*/

	/* Problem specific variables（问题具体变量）*/
	//要优化的问题的参数数量
	int D = 100; /*The number of parameters of the problem to be optimized*/
	double lb = -5.12; /*lower bound of the parameters. */
	double ub = 5.12; /*upper bound of the parameters. lb and ub can be defined as arrays for the problems of which parameters have different bounds（lb和ub可以被定义为参数具有不同边界的问题的数组）*/

	int runtime = 30;  /*Algorithm can be run many times in order to see its robustness（算法可以运行多次以查看其鲁棒性）*/

	//int dizi1[]=new int[10];
	double Foods[][]=new double[FoodNumber][D];        /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber（食物是食物源的人口。 食物矩阵的每一行都是一个保存要优化的D参数的向量。 食物矩阵的行数等于FoodNumber）*/
	//就是Y值
	double f[]=new double[FoodNumber];        /*f is a vector holding objective function values associated with food sources（f是持有与食物来源相关的目标函数值的向量） */
	double fitness[]=new double[FoodNumber];      /*fitness is a vector holding fitness (quality) values associated with food sources（fitness是持有与食物来源相关的fitness（质量）值的向量）*/
	double trial[]=new double[FoodNumber];         /*trial is a vector holding trial numbers through which solutions can not be improved（trial是一个持有trial数大小并解决方案在trial次数之后无法改进的向量）*/
	double prob[]=new double[FoodNumber];          /*prob is a vector holding probabilities of food sources (solutions) to be chosen（prob是一个食物源（解决方案）被选择的概率的向量）*/
	//就是X值
	double solution[]=new double[D];            /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
	
	
    //新解决方案的目标函数值               
	double ObjValSol;              /*Objective function value of new solution*/
	//新解决方案的fitness值
	double FitnessSol;              /*Fitness value of new solution*/
	int neighbour, param2change;                   /*param2change corrresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/

	//ABC算法得到的最优解
	double GlobalMin;                       /*Optimum solution obtained by ABC algorithm*/
	//最优解的参数
	double GlobalParams[]=new double[D];                   /*Parameters of the optimum solution*/
	//GlobalMins存放多次运行中每次运行的GlobalMin
	double GlobalMins[]=new double[runtime];            
	         /*GlobalMins holds the GlobalMin of each run in multiple runs*/
	//范围在[0,1)之间的随机数
	double r; /*a random number in the range [0,1)*/

	/*a function pointer returning double and taking a D-dimensional array as argument */
	/*If your function takes additional arguments then change function pointer definition and lines calling "...=function(solution);" in the code*/


//	typedef double (*FunctionCallback)(double sol[D]);  

	/*benchmark functions */

//	double sphere(double sol[D]);
//	double Rosenbrock(double sol[D]);
//	double Griewank(double sol[D]);
//	double Rastrigin(double sol[D]);

	/*Write your own objective function name instead of sphere*/
//	FunctionCallback function = &sphere;

	/*Fitness function*/
	double CalculateFitness(double fun) 
	{
		 double result=0;
		 if(fun>=0)
		 {
			 result=1/(fun+1);
		 }
		 else
		 {
			 result=1+Math.abs(fun);
		 }
		 return result;
	}

	/*
	 * The best food source is memorized
	 * 记忆最好的食物来源
	 * */
	void MemorizeBestSource() 
	{
		int i,j;
	
		for(i=0;i<FoodNumber;i++)
		{
			if (f[i]<GlobalMin)
			{
				GlobalMin=f[i];
				for(j=0;j<D;j++)
				GlobalParams[j]=Foods[i][j];
			}
		}
	}

	/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */
	/* Counters of food sources are also initialized in this function*/
	
	void init(int index)
	{
	    int j;
	    for (j=0;j<D;j++)
	    {
	    	//Math.random()返回的是double范围的[0.0,1.0)
	    	//2^15=32768
			r = (   (double)Math.random()*32767 / ((double)32767+(double)(1)) );
			Foods[index][j]=r*(ub-lb)+lb;
			//solution[j]其实就是可行解的值，x值,也就是说FoodNumber越大，可行解的值越多，区间范围寻找越广.这个想法还不知对不对，应该是对的
			solution[j]=Foods[index][j];
		}
		f[index]=calculateFunction(solution);
		fitness[index]=CalculateFitness(f[index]);
		trial[index]=0;
	}


	/*All food sources are initialized */
	void initial()
	{
		int i;
		for(i=0;i<FoodNumber;i++)
		{
			init(i);
		}
		GlobalMin=f[0];
	    for(i=0;i<D;i++)
	    GlobalParams[i]=Foods[0][i];
	}

	//派遣采蜜蜂
	void SendEmployedBees()
	{
	    int i,j;
	    /*Employed Bee Phase*/
	    for (i=0;i<FoodNumber;i++)
		{
			/*The parameter to be changed is determined randomly*/
			r = ((double) Math.random()*32767 / ((double)(32767)+(double)(1)) );
			param2change=(int)(r*D);
			
			/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
			r = ((double)Math.random()*32767 / ((double)(32767)+(double)(1)) );
			neighbour=(int)(r*FoodNumber);
			
			/*Randomly selected solution must be different from the solution i*/        
		    // while(neighbour==i)
		    // {
		    // r = (   (double)Math.random()*32767 / ((double)(32767)+(double)(1)) );
		    // neighbour=(int)(r*FoodNumber);
		    // }
			for(j=0;j<D;j++)
			solution[j]=Foods[i][j];

			/*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
			r =((double)Math.random()*32767 / ((double)(32767)+(double)(1)));
			solution[param2change]=Foods[i][param2change]+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2;

			/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
			if (solution[param2change]<lb)
			   solution[param2change]=lb;
			if (solution[param2change]>ub)
			   solution[param2change]=ub;
			ObjValSol=calculateFunction(solution);
			FitnessSol=CalculateFitness(ObjValSol);
			
			/*	a greedy selection is applied between the current solution i and its mutant
			 * 	在当前的解决方案i和它的突变体之间应用贪婪的选择
			 */
			if (FitnessSol>fitness[i])
			{
				/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
				trial[i]=0;
				for(j=0;j<D;j++)
				Foods[i][j]=solution[j];
				f[i]=ObjValSol;
				fitness[i]=FitnessSol;
			}
			else
			{
				/*if the solution i can not be improved, increase its trial counter*/
				trial[i]=trial[i]+1;
			}
			
		}
		/*end of employed bee phase*/
	}

	/* A food source is chosen with the probability which is proportioal to its quality*/
	/*Different schemes can be used to calculate the probability values*/
	/*For example prob(i)=fitness(i)/sum(fitness)*/
	/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
	/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/
	void CalculateProbabilities()
	{
	     int i;
	     double maxfit;
	     maxfit=fitness[0];
	     for (i=1;i<FoodNumber;i++)
	     {
	         if (fitness[i]>maxfit)
	         maxfit=fitness[i];
	     }

	     for (i=0;i<FoodNumber;i++)
	     {
	    	 prob[i]=(0.9*(fitness[i]/maxfit))+0.1;
	     }

	}

	//派遣观察蜂
	void SendOnlookerBees()
	{
	  int i,j,t;
	  i=0;
	  t=0;
	  /*onlooker Bee Phase*/
	  while(t<FoodNumber)
	  {
        r = ((double)Math.random()*32767 / ((double)(32767)+(double)(1)) );
        if(r<prob[i]) /*choose a food source depending on its probability to be chosen*/
        {
			t++;
			/*The parameter to be changed is determined randomly*/
			r = ((double)Math.random()*32767 / ((double)(32767)+(double)(1)) );
			param2change=(int)(r*D);
			
			/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
			r = ((double)Math.random()*32767 / ((double)(32767)+(double)(1)) );
			neighbour=(int)(r*FoodNumber);

			/*
			 * Randomly selected solution must be different from the solution i
			 * 随机选择的解决方案必须与解决方案i不同
			 * */        
			while(neighbour == i)
			{
				//System.out.println(Math.random()*32767+"  "+32767);
				r = ((double)Math.random()*32767 / ((double)(32767)+(double)(1)) );
				neighbour=(int)(r*FoodNumber);
			}
			for(j=0;j<D;j++)
			solution[j]=Foods[i][j];

			/*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
			r = (   (double)Math.random()*32767 / ((double)(32767)+(double)(1)) );
			solution[param2change]=Foods[i][param2change]+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2;

			/*
			 * if generated parameter value is out of boundaries, it is shifted onto the boundaries
			 * 如果生成的参数值超出边界，则会移至边界
			 */
			if (solution[param2change]<lb)
			   solution[param2change]=lb;
			if (solution[param2change]>ub)
			   solution[param2change]=ub;
			ObjValSol=calculateFunction(solution);
			FitnessSol=CalculateFitness(ObjValSol);
			
			/*
			 * a greedy selection is applied between the current solution i and its mutant
			 * 在当前的解决方案i和它的突变体之间应用贪婪的选择
			 * */
			if (FitnessSol>fitness[i])
			{
				/*
				 * If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i
				 * 如果突变体解决方案比当前的解决方案更好，则用突变体替换解决方案，并重置解决方案i的试验计数器
				 * */
				trial[i]=0;
				for(j=0;j<D;j++)
				Foods[i][j]=solution[j];
				f[i]=ObjValSol;
				fitness[i]=FitnessSol;
			}
			else
			{   
				/*
				 * if the solution i can not be improved, increase its trial counter
				 * 如果解决方案i无法改进，请增加其试用计数器
				 * */
				trial[i]=trial[i]+1;
			}
        }/*if */
        i++;
        if (i==FoodNumber)
        i=0;
	  }/*while*/
	  /*end of onlooker bee phase     */
	}

	/*
	 * determine the food sources whose trial counter exceeds the "limit" value. In Basic ABC, only one scout is allowed to occur in each cycle
	 * 确定试验计数器超过“极限”值的食物来源。 在基本ABC中，每个周期只允许发生一次侦察
	 * */
	//派遣侦查蜂
	void SendScoutBees()
	{
		int maxtrialindex,i;
		maxtrialindex=0;
		for (i=1;i<FoodNumber;i++)
		{
			if (trial[i]>trial[maxtrialindex])
			maxtrialindex=i;
		}
		if(trial[maxtrialindex]>=limit)
		{
			init(maxtrialindex);
		}
	}

	double calculateFunction(double sol[])
	{
		return Rastrigin (sol);	
	}
	
	//领域
	double sphere(double sol[])
	{
		int j;
		double top=0;
		for(j=0;j<D;j++)
		{
			top=top+sol[j]*sol[j];
		}
		return top;
	}

	double Rosenbrock(double sol[])
	{
		int j;
		double top=0;
		for(j=0;j<D-1;j++)
		{
			top=top+100*Math.pow((sol[j+1]-Math.pow((sol[j]),(double)2)),(double)2)+Math.pow((sol[j]-1),(double)2);
		}
		return top;
	}

	 double Griewank(double sol[])
	 {
		 int j;
		 double top1,top2,top;
		 top=0;
		 top1=0;
		 top2=1;
		 for(j=0;j<D;j++)
		 {
			top1=top1+Math.pow((sol[j]),(double)2);
			top2=top2*Math.cos((((sol[j])/Math.sqrt((double)(j+1)))*Math.PI)/180);
		 }	
		 top=(1/(double)4000)*top1-top2+1;
		 return top;
	 }

	 double Rastrigin(double sol[])
	 {
		 int j;
		 double top=0;
		 //Math.pow(底数,几次方) http://blog.csdn.net/zhview/article/details/53198581
		 //cos() 方法用于返回指定double类型参数的余弦值。
		 for(j=0;j<D;j++)
		 {
			 //求这个函数的极大值（最优解）
			top=top+(Math.pow(sol[j],(double)2)-10*Math.cos(2*Math.PI*sol[j])+10);
		 }
		 return top;
	 }
}
