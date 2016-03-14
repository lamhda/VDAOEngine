package vdaoengine.analysis.grammars;

import java.io.File;
import java.util.Vector;

import vdaoengine.analysis.elmap.ElmapAlgorithmEpoch;
import vdaoengine.data.VDataSet;
import vdaoengine.data.VDataTable;
import vdaoengine.data.io.VDatReadWrite;
import vdaoengine.utils.OptionParser;
import vdaoengine.utils.VSimpleProcedures;

public class ComputePrincipalGraph{

	public ConfigFile config = new ConfigFile();
	public VDataSet dataset = null;
	public String project = "";

	public ElasticEnergyOptimization elo = null;
	
	public BaseOptimizationAlgorithm alg = null;

	
	public Graph graph = null;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try{
			
			if(args.length==0){
				args = new String[8];
				args[0] = "--dat";
				//args[1] = "c:/datas/elastictree/tree23.dat";
				//args[1] = "c:/datas/elastictree/iris.dat";
				//args[1] = "c:/datas/elastictree/wine.dat";
				//args[1] = "c:/datas/elastictree/wdbc_cut.dat";
				//args[1] = "c:/datas/elastictree/winequality-red.dat";
				//args[1] = "c:/datas/elastictree/scoelicolor.dat";
				//args[1] = "c:/datas/elastictree/forestfires.dat";
				//args[1] = "c:/datas/elastictree/abalone.dat";
				args[1] = "c:/datas/elastictree/wangn5000_t.dat";
				args[2] = "--config";
				args[3] = "elmap.ini";
				args[4] = "--num";
				args[5] = "14";
				args[6] = "--normalize";
				args[7] = "1";
			}
			
			OptionParser options = new OptionParser(args, null);
			String s = null;
			File f = null;
			File datFile = null;
			File configFile = null;
			int algorithmNumber = 0;
			boolean normalize = false;
			Boolean b = false;
			
			if ((f = options.fileRequiredOption("dat", "dat file")) != null)
				datFile = f;
			if ((f = options.fileRequiredOption("config", "configuration file")) != null)
				configFile = f;
			if ((s = options.stringOption("num", "number of the algorithm in the config file")) != null)
				algorithmNumber = Integer.parseInt(s);
			if ((s = options.stringOption("normalize", "do standard data normalization")) != null)
				normalize = s.equals("1");
			options.done();
			
			ComputePrincipalGraph cpg = new ComputePrincipalGraph();
			String projectName = datFile.getName().substring(0,datFile.getName().length()-4);
			System.out.println("Project "+projectName);
			cpg.project = projectName;
			cpg.config.readFile(configFile.getAbsolutePath(), algorithmNumber);
			
			VDataTable vt = VDatReadWrite.LoadFromVDatFile(datFile.getAbsolutePath());
			if(!normalize) 
				cpg.dataset = VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(vt, -1);
			else
				cpg.dataset = VSimpleProcedures.SimplyPreparedDataset(vt, -1);
			
			cpg.compute();
			
			
			cpg.saveToFile(datFile.getAbsolutePath().substring(0,datFile.getAbsolutePath().length()-4)+".vem");
						
			cpg.graph.writeOutNodes(datFile.getAbsolutePath().substring(0,datFile.getAbsolutePath().length()-4)+".nodes");
			cpg.graph.writeOutEdges(datFile.getAbsolutePath().substring(0,datFile.getAbsolutePath().length()-4)+".edges");
			
			vdaoengine.data.io.VDatReadWrite.useQuotesEverywhere = false;
			vdaoengine.data.io.VDatReadWrite.saveToSimpleDatFilePureNumerical(cpg.dataset,datFile.getAbsolutePath().substring(0,datFile.getAbsolutePath().length()-4)+".data");
			
			
		}catch(Exception e){
			e.printStackTrace();
		}

	}
	
	public void compute() throws Exception{
		
		dataset.calcStatistics();
		System.out.println("Variation = "+dataset.simpleStatistics.totalDispersion);
		dataset.simpleStatistics.calculate();
		System.out.println("Variation = "+dataset.simpleStatistics.totalDispersion);
		
		
		graph = new Graph();
		
		elo = new ElasticEnergyOptimization(dataset, graph);
		alg = new BaseOptimizationAlgorithm(dataset);
		alg.verbose = false;
		
		//alg.setElementFeeExponential(1e-4f, 3f);
		//alg.setElementFeeLinear(0.3e-4f, 1f);

		defineGrammarType();
		
		for(int i=0;i<config.epochs.size();i++){
			System.out.println("Epoch: "+(i+1));
			ElmapAlgorithmEpoch ep = config.epochs.get(i);
			graph.setDefaultEdgeElasticityCoeff(ep.EP);
			graph.setDefaultElasticityCoeffs(ep.RP);		
			
			alg.parameters.maxNumberOfNodes = ep.numberOfIterations;
			//alg.movieFolder = "c:/datas/elastictree/movie/";
			alg.setGraph(graph);
			alg.initializeGraph();
			
			if(ep.minimize){
				alg.run(elo);
				if(alg.convergedByComplexity)
					System.out.println("Converged by complexity");
				
				elo.updateEnergyValue();
				
				System.out.println("MSE = "+alg.graph.calcMSE(dataset, elo.taxons));
				System.out.println("Energy = "+elo.energyValue);
			}
		
		}
		graph = alg.graph;
	}
	
	public void defineGrammarType() throws Exception{
		
		alg.grammars = new Vector<GraphGrammar>();
		
		if(config.grammartype.equals("tree")){
			GraphGrammar grammarGrow = new GraphGrammar(); 
			BisectEdge be = new BisectEdge();
			AddNodeToNode an = new AddNodeToNode(); 
			grammarGrow.operations.add(be);
			grammarGrow.operations.add(an);		
			alg.grammars.add(grammarGrow);
		}
		if(config.grammartype.equals("treeWithTrimming")){
			GraphGrammar grammarGrow = new GraphGrammar(); 
			BisectEdge be = new BisectEdge();
			AddNodeToNode an = new AddNodeToNode(); 
			grammarGrow.operations.add(be);
			grammarGrow.operations.add(an);			
			GraphGrammar grammarShrink = new GraphGrammar();
			RemoveLeaf rl = new RemoveLeaf();
			RemoveInternalEdge rie = new RemoveInternalEdge(); 
			grammarShrink.operations.add(rl);
			grammarShrink.operations.add(rie);
			alg.grammars.add(grammarGrow);
			alg.grammars.add(grammarGrow);
			alg.grammars.add(grammarShrink);
		}
	}
	
	public void saveToFile(String fn) throws Exception{
		alg.graph.recalcIndexMaps();
		
		alg.graph.saveToFile(fn,project);
		
	}

}
