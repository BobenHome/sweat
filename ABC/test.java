package test;

public class test
{
    
    static beeColony bee = new beeColony();
    
    public static void main(String[] args)
    {
        int iter = 0;
        int run = 0;
        int j = 0;
        double mean = 0;
        //srand(time(NULL));
        for (run = 0; run < bee.runtime; run++)
        {
            bee.initial();
            bee.MemorizeBestSource();
            for (iter = 0; iter < bee.maxCycle; iter++)
            {
                bee.SendEmployedBees();
                bee.CalculateProbabilities();
                bee.SendOnlookerBees();
                bee.MemorizeBestSource();
                bee.SendScoutBees();
            }
            for (j = 0; j < bee.D; j++)
            {
                //System.out.println("GlobalParam[%d]: %f\n",j+1,GlobalParams[j]);
                System.out.println("GlobalParam[" + (j + 1) + "]:" + bee.GlobalParams[j]);
            }
            //System.out.println("%d. run: %e \n",run+1,GlobalMin);
            System.out.println((run + 1) + ".run:" + bee.GlobalMin);
            bee.GlobalMins[run] = bee.GlobalMin;
            mean = mean + bee.GlobalMin;
        }
        mean = mean / bee.runtime;
        //System.out.println("Means of %d runs: %e\n",runtime,mean);
        System.out.println("Means of " + bee.runtime + "runs: " + mean);
        //System.out.println("极大值："+(Math.pow(mean,(double)2)-10*Math.cos(2*Math.PI*mean)+10));
        for (int i = 0; i < bee.GlobalMins.length - 1; i++)
        {
            for (int k = 0; k < bee.GlobalMins.length - 1 - i; k++)
            {
                if (bee.GlobalMins[k] > bee.GlobalMins[k + 1])
                {
                    double temp = bee.GlobalMins[k];
                    bee.GlobalMins[k] = bee.GlobalMins[k + 1];
                    bee.GlobalMins[k + 1] = temp;
                }
            }
        }
        
        for (int i = 0; i < bee.GlobalMins.length; i++)
        {
            System.out.println(i + 1 + " " + bee.GlobalMins[i]);
        }
        
    }
    
}
