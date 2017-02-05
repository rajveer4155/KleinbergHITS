import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class Hits {
	
	public static String file;
	public static int iterations=-1;
	public static double initialVal;
	public static int numOfVertices;
	public static int numOfEdges;
	public static double[] a; // authority values
	public static double[] h; // hub values
	public static double[] preva; // Previous authority values
	public static double[] prevh; // Previous hub values
	public double errorRate=0.00005;
	DirectedGraph<Integer,Integer> dgraph;
	
	public Hits(){}
	
	public void initializeGraph() throws FileNotFoundException, IOException{
		
		//** Reading Graph Data from file for initialization

		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
		    String line;
		    int i=0;
		   
		    while ((line = br.readLine()) != null) {
		       // process the line.
		    	String[] str=line.split(" ");
		    	
		    	if(i==0){ //If first line
		    		numOfVertices=Integer.parseInt(str[0]);
		    		System.out.println("Num of Vertices Read: "+numOfVertices);
		    		numOfEdges=Integer.parseInt(str[1]);
		    		System.out.println("Num of Edges Read: "+numOfEdges);
		    		//Also Initializing AdjMatrix
		    		dgraph=new DirectedGraph<Integer,Integer>(numOfVertices);
		    		dgraph.initializeAdjMatrices(numOfVertices);
		    		
		    	}else{ // For all other Lines
		    		
		    		int u=Integer.parseInt(str[0]);
		    		int v=Integer.parseInt(str[1]);
		    		System.out.println(u+"-"+v);
		    		
		    		//** Adding Vertices and Edge in Adj List
		    		dgraph.addEdge(dgraph.addVertice(u, null),dgraph.addVertice(v, null),-1);
		    		
		    		//** Adding Edge in AdjMatrix
		    		dgraph.setEdgeInAdjMatrix(u, v);
		    	}
		    	i++;
		    } //while ends
		} // try ends	
		System.out.println("Num of Vertices: "+dgraph.ESet.size());
		
		System.out.println("Num of Edges : "+dgraph.VSet.size());
		
	}
	public void initializeKleinberg(){
		//Initializing Hub and Authority Arrays
		a = new double[numOfVertices];
		h = new double[numOfVertices];
		
		preva = new double[numOfVertices];
		prevh = new double[numOfVertices];
		
		if(numOfVertices>10){
			initialVal=-1.00000;
		}
		
		if(initialVal==0){
			initialVal=0.00000;
		}else if(initialVal==1){
			initialVal=1.00000;
		}else if(initialVal==-1){
			initialVal=(1.00000)/numOfVertices;
		}else if(initialVal==-2){
			initialVal=(1.00000)/(Math.sqrt(numOfVertices));
		}
		
		//Initializing Authority and Hub Values
		//double scaleFactor=Math.sqrt((numOfVertices*initialVal));
		for(int i=0;i<numOfVertices;i++){
			//a[i]=initialVal/scaleFactor;
			//h[i]=initialVal/scaleFactor;
			a[i]=initialVal;
			h[i]=initialVal;
		}
		
	}
	
	//*** Checking Convergence
	public boolean checkConvergence(){
		boolean flag=true;
		for(int i=0;i<numOfVertices;i++){
			if((preva[i]-a[i])>errorRate || (prevh[i]-h[i])>errorRate){
				flag=false;
			}
		} 
		return flag;
	}
	
	//*** Starting KleinBerg
	public void runKleinbergHits(){
		int iter=0;
		System.out.print("Base : 0 :");
		for(int i = 0; i < numOfVertices; i++) {
            System.out.print(" A/H[ "+i+"]="+a[i]+"/"+h[i]); 
        }
		System.out.println();
		do{
			++iter;
			System.out.println("Iteration: "+iter);
			
			//Assigning Prev A and Prev H vals
			for(int i=0;i<numOfVertices;i++){
				preva[i]=a[i];
				prevh[i]=h[i];
			}
			
			//Perform authority update
			double scaleFactor=0;
			for(int i=0;i<numOfVertices;i++){
				
				for(int j=0;j<numOfVertices;j++){
					if(dgraph.isEdgeInAdjMatrix(j,i)==1){ //checking transpose edge
						a[i]+=prevh[j];
					}
				}
				scaleFactor+=(a[i]*a[i]);
			} // Authority update ends
			scaleFactor=Math.sqrt(scaleFactor);
			
			//Performing scaling on A-Update
			for(int i=0;i<numOfVertices;i++){
				a[i]=a[i]/scaleFactor;
			} 

			
			//Perform Hub-update
		    scaleFactor=0;
			for(int i=0;i<numOfVertices;i++){
				for(int j=0;j<numOfVertices;j++){
					if(dgraph.isEdgeInAdjMatrix(j,i)==1){ //checking transpose edge
						h[i]+=a[j];
					}
				}
				scaleFactor+=(h[i]*h[i]);
			} //Hub-update ends
			
			scaleFactor=Math.sqrt(scaleFactor);
			//Performing scaling on A-Update
			for(int i=0;i<numOfVertices;i++){
				h[i]=h[i]/scaleFactor;
			} 
			
			System.out.print("Iter : "+iter+" :");
			for(int i = 0; i < numOfVertices; i++) {
	            System.out.print(" A/H[ "+i+"]="+a[i]+"/"+h[i]); 
	        }
			System.out.println();
			
		}while(!checkConvergence());
		
		System.out.println("Converged on Iteration: "+iter);
	}
	
	
	//*** Main Function
	public static void main(String args[]) throws IOException{
		
		  if (args.length < 1) {
	            System.out.println("Error!Please provide file name.");
	            System.out.println("Example to run:");
	            System.out.println("hits <iterations> <initialvalue> <graphfilename>");
	            return;
	        }
	        if (!new File(args[2]).exists()) {
	            System.out.println("File not found. Please enter an existing file name.");
	            return;
	        }
	     
	     //*** Reading Values
	        
	     iterations=Integer.parseInt(args[0]);
	 	 initialVal=Double.parseDouble(args[1]);
	 	 file=args[2];
	 	 
	     //*** Checking if file is empty
	     File fr = new File(file);
	     if(fr.length()==0){
	    	 System.out.println("Error! File is Empty");
	    	 return;
	     }
	     
	     Hits hits=new Hits();
	     hits.initializeGraph();
	     hits.initializeKleinberg();
	     hits.runKleinbergHits();
			
	} // Main ends
}
