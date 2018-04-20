package beeColony;

public class test {
	
	static beeColony bee=new beeColony();
	public static void main(String[] args) {
		int iter=0;
		int run=0;
		int j=0;
		double mean=0;
		//srand(time(NULL));
		for(run=0;run<bee.runtime;run++)
		{
			bee.initial();
			bee.MemorizeBestSource();
			for (iter=0;iter<bee.maxCycle;iter++)
			{
				bee.SendEmployedBees();
				bee.CalculateProbabilities();
				bee.SendOnlookerBees();
				bee.MemorizeBestSource();
				bee.SendScoutBees();
			}
			for(j=0;j<bee.D;j++)
			{
				//System.out.println("GlobalParam[%d]: %f\n",j+1,GlobalParams[j]);
				System.out.println("GlobalParam["+(j+1)+"]:"+bee.GlobalParams[j]);
			}
			//System.out.println("%d. run: %e \n",run+1,GlobalMin);
			System.out.println((run+1)+".run:"+bee.GlobalMin);
			bee.GlobalMins[run]=bee.GlobalMin;
			mean=mean+bee.GlobalMin;
		}
		mean=mean/bee.runtime;
		//System.out.println("Means of %d runs: %e\n",runtime,mean);
		System.out.println("Means of "+bee.runtime+"runs: "+mean);
		System.out.println("极大值："+(Math.pow(mean,(double)2)-10*Math.cos(2*Math.PI*mean)+10));
		
	}

}
