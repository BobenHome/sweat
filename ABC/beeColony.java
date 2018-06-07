package test;

import java.lang.Math;

public class beeColony
{
    /*ABC算法的控制参数*/
    //蜂群大小（采蜜蜂数量+观察蜂数量）
    int NP = 20;
    
    //食物源大小（可行解的数量）
    int FoodNumber = NP / 2;
    
    //无法通过limit次试验来改进的食物源将被采蜜蜂遗弃
    int limit = 100;
    
    //觅食的周期数{停止标准}
    int maxCycle = 2500;
    
    /*问题具体变量*/
    //要优化的问题的参数数量
    int D = 100;
    
    /*lb和ub可以被定义为参数具有不同边界的问题的数组*/
    double lb = -5.12;
    
    double ub = 5.12;
    
    /*算法可以运行多次以查看其鲁棒性*/
    int runtime = 30;
    
    /*食物是食物源的人口。 食物矩阵的每一行都是一个保存要优化的D参数的向量。 食物矩阵的行数等于FoodNumber*/
    double Foods[][] = new double[FoodNumber][D];
    
    /*f是持有与食物来源相关的目标函数值的向量 */
    //就是Y值
    double f[] = new double[FoodNumber];
    
    /*fitness是持有与食物来源相关的fitness（质量）值的向量*/
    double fitness[] = new double[FoodNumber];
    
    /*trial是一个持有trial数大小并解决方案在trial次数之后无法改进的向量*/
    double trial[] = new double[FoodNumber];
    
    /*prob是一个食物源（解决方案）被选择的概率的向量*/
    double prob[] = new double[FoodNumber];
    
    //就是X值
    double solution[] = new double[D];
    
    //新解决方案的目标函数值               
    double ObjValSol; /*Objective function value of new solution*/
    
    //新解决方案的fitness值
    double FitnessSol;
    
    int neighbour, param2change;
    
    //ABC算法得到的最优解
    double GlobalMin;
    
    //最优解的参数
    double GlobalParams[] = new double[D];
    
    //GlobalMins存放多次运行中每次运行的GlobalMin
    double GlobalMins[] = new double[runtime];
    
    //范围在[0,1)之间的随机数
    double r;
    
    //与全局最优值交叉时搜索领域系数
    double cr = 0.6;
    
    int bestLocationC = 3;
    
    /*基准函数*/
    //double sphere(double sol[D]);
    //double Rosenbrock(double sol[D]);
    //double Griewank(double sol[D]);
    //double Rastrigin(double sol[D]);
    
    /*计算fitness值*/
    double CalculateFitness(double fun)
    {
        double result = 0;
        if (fun >= 0)
        {
            result = 1 / (fun + 1);
        }
        else
        {
            result = 1 + Math.abs(fun);
        }
        return result;
    }
    
    /*
     * 记忆最好的食物来源
     * 通过比较f[]的值
     * */
    void MemorizeBestSource()
    {
        int i, j;
        
        for (i = 0; i < FoodNumber; i++)
        {
            //取得函数值最小
            if (f[i] < GlobalMin)
            {
                GlobalMin = f[i];
                for (j = 0; j < D; j++)
                    GlobalParams[j] = Foods[i][j];
            }
        }
    }
    
    /**
     * 
     * <功能详细描述>
     * @param index:可行解序号
     * 初始化
     */
    void init(int index)
    {
        int j;
        for (j = 0; j < D; j++)
        {
            //Math.random()返回的是double范围的[0.0,1.0)
            //2^15=32768
            r = ((double)Math.random() * 32767 / ((double)32767 + (double)(1)));
            //Foods[][]存放的值区间范围都是[-5.12,5.12]
            Foods[index][j] = r * (ub - lb) + lb;
            //solution[j]其实就是可行解的值，x值,也就是说FoodNumber越大，可行解的值越多，区间范围寻找越广.这个想法还不知对不对，应该是对的
            //上面这个想法不对，区间范围还是在[-5.12,5.12]，当D的值越大，表明当一个可行解的时候（也就是index），测试的值越多。
            solution[j] = Foods[index][j];
        }
        //f[]存放的就是每个可行解的假Y值（多个D参数测试值累计加起来的值）
        f[index] = calculateFunction(solution);
        //fitness[]存放的健康值，但是健康值的计算却让我一头雾水
        //fitness越小表示函数值越大，也就是极大值，fitness越大表示函数值越小，也就是极小值
        fitness[index] = CalculateFitness(f[index]);
        trial[index] = 0;
    }
    
    /**
     * <功能详细描述>
     * 所有 food sources 被初始化
     */
    void initial()
    {
        int i;
        for (i = 0; i < FoodNumber; i++)
        {
            init(i);
        }
        //GlobalMin全局最优值，初始化全局最优值，假设第一个数为最优
        GlobalMin = f[0];
        //GlobalParams[]存放全局最优参数，假设第一个可行解的参数为最优
        for (i = 0; i < D; i++)
            GlobalParams[i] = Foods[0][i];
    }
    
    //派遣采蜜蜂
    void SendEmployedBees()
    {
        int i, j;
        /*采蜜蜂开始*/
        for (i = 0; i < FoodNumber; i++)
        {
            //要更改的参数是随机确定的
            r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
            //param2change为[0,100)之间的值
            param2change = (int)(r * D);
            
            //随机选择的solution是相对于solution[i]的突变solution
            r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
            //neighbour为[0,10)之间的值
            neighbour = (int)(r * FoodNumber);
            
            //随机选择的解决方案必须与解决方案i不同
            while (neighbour == i)
            {
                r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
                neighbour = (int)(r * FoodNumber);
            }
            for (j = 0; j < D; j++)
                solution[j] = Foods[i][j];
            
            r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
            //(r-0.5)*2 范围在[-1,1)之间
            //[-5.12,5.12]-[-5.12,5.12]范围在[-10.24,10.24]之间
            solution[param2change] =
                Foods[i][param2change] + (Foods[i][param2change] - Foods[neighbour][param2change]) * (r - 0.5) * 2;
            
            
            
            //如果生成的参数值超出边界，则会移至边界
            if (solution[param2change] < lb)
                solution[param2change] = lb;
            if (solution[param2change] > ub)
                solution[param2change] = ub;
            ObjValSol = calculateFunction(solution);
            FitnessSol = CalculateFitness(ObjValSol);
            
            //在当前的解决方案i和它的突变体之间应用贪婪的选择
            if (FitnessSol > fitness[i])
            {
                //如果突变体解决方案比当前的解决方案更好，则用突变体替换解决方案，并重置解决方案i的试验计数器
                trial[i] = 0;
                for (j = 0; j < D; j++)
                    Foods[i][j] = solution[j];
                f[i] = ObjValSol;
                fitness[i] = FitnessSol;
            }
            else
            {
                //如果解决方案i无法改进，请增加其试用计数器
                trial[i] = trial[i] + 1;
            }
            
        }
        /*采蜜锋结束*/
    }
    
    /*选择食物来源的概率与其质量成正比*/
    /*可以使用不同的方案来计算概率值*/
    /*例如prob(i)=fitness(i)/sum(fitness)*/
    /*或者使用这种方式 prob(i)=a*fitness(i)/max(fitness)+b*/
    /*通过使用适应值来计算概率值，并通过除以最大适应值来归一化*/
    void CalculateProbabilities()
    {
        int i;
        double maxfit;
        maxfit = fitness[0];
        for (i = 1; i < FoodNumber; i++)
        {
            if (fitness[i] > maxfit)
                maxfit = fitness[i];
        }
        
        for (i = 0; i < FoodNumber; i++)
        {
            prob[i] = (0.9 * (fitness[i] / maxfit)) + 0.1;
        }
        
    }
    
    //派遣观察蜂
    void SendOnlookerBees()
    {
        int i, j, t;
        i = 0;
        t = 0;
        /*观察蜂开始*/
        while (t < FoodNumber)
        {
            r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
            //根据选择的概率选择食物来源0.68
            if (0.68 < prob[i])
            {
                t++;
                
                r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
                param2change = (int)(r * D);
                
                r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
                neighbour = (int)(r * FoodNumber);
                
                //随机选择的解决方案必须与解决方案i不同
                while (neighbour == i)
                {
                    r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
                    neighbour = (int)(r * FoodNumber);
                }
                for (j = 0; j < D; j++)
                    solution[j] = Foods[i][j];
                
                r = ((double)Math.random() * 32767 / ((double)(32767) + (double)(1)));
                solution[param2change] =
                    Foods[i][param2change] + (Foods[i][param2change] - Foods[neighbour][param2change]) * (r - 0.5) * 2;
                
                //如果生成的参数值超出边界，则会移至边界
                if (solution[param2change] < lb)
                    solution[param2change] = lb;
                if (solution[param2change] > ub)
                    solution[param2change] = ub;
                ObjValSol = calculateFunction(solution);
                FitnessSol = CalculateFitness(ObjValSol);
                
                //在当前的解决方案i和它的突变体之间应用贪婪的选择
                if (FitnessSol > fitness[i])
                {
                    //如果突变体解决方案比当前的解决方案更好，则用突变体替换解决方案，并重置解决方案i的试验计数器
                    trial[i] = 0;
                    for (j = 0; j < D; j++)
                        Foods[i][j] = solution[j];
                    f[i] = ObjValSol;
                    fitness[i] = FitnessSol;
                }
                else
                {
                    //如果解决方案i无法改进，请增加其试用计数器
                    trial[i] = trial[i] + 1;
                }
            }
            i++;
            if (i == FoodNumber)
                i = 0;
        }
        /*观察蜂结束*/
    }
    
    //确定试验计数器超过"极限值"的食物来源。 在基本ABC中，每个周期只允许发生一次侦察
    //派遣侦查蜂
    void SendScoutBees()
    {
        int maxtrialindex, i;
        maxtrialindex = 0;
        for (i = 1; i < FoodNumber; i++)
        {
            if (trial[i] > trial[maxtrialindex])
                maxtrialindex = i;
        }
        if (trial[maxtrialindex] >= limit)
        {
            init(maxtrialindex);
        }
    }
    
    double calculateFunction(double sol[])
    {
        return Rastrigin(sol);
    }
    
    double sphere(double sol[])
    {
        int j;
        double top = 0;
        for (j = 0; j < D; j++)
        {
            top = top + sol[j] * sol[j];
        }
        return top;
    }
    
    double Rosenbrock(double sol[])
    {
        int j;
        double top = 0;
        for (j = 0; j < D - 1; j++)
        {
            top = top + 100 * Math.pow((sol[j + 1] - Math.pow((sol[j]), (double)2)), (double)2)
                + Math.pow((sol[j] - 1), (double)2);
        }
        return top;
    }
    
    double Griewank(double sol[])
    {
        int j;
        double top1, top2, top;
        top = 0;
        top1 = 0;
        top2 = 1;
        for (j = 0; j < D; j++)
        {
            top1 = top1 + Math.pow((sol[j]), (double)2);
            top2 = top2 * Math.cos((((sol[j]) / Math.sqrt((double)(j + 1))) * Math.PI) / 180);
        }
        top = (1 / (double)4000) * top1 - top2 + 1;
        return top;
    }
    
    //该函数计算了一个可行解时多个D参数测试值累计的Y值
    double Rastrigin(double sol[])
    {
        int j;
        double top = 0;
        //Math.pow(底数,几次方) http://blog.csdn.net/zhview/article/details/53198581
        //Math.cos() 方法用于返回指定double类型参数的余弦值。
        for (j = 0; j < D; j++)
        {
            //求这个函数的极值（最优解）
            top = top + (Math.pow(sol[j], (double)2) - 10 * Math.cos(2 * Math.PI * sol[j]) + 10);
        }
        return top;
    }
}
