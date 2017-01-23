package vdaoengine;

import java.io.File;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.HashMap;
import java.util.Locale;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;

import vdaoengine.analysis.PCAMethod;
import vdaoengine.data.VDataSet;
import vdaoengine.data.VDataTable;
import vdaoengine.data.io.VDatReadWrite;
import vdaoengine.utils.Algorithms;
import vdaoengine.utils.MetaGene;
import vdaoengine.utils.Utils;
import vdaoengine.utils.VSimpleFunctions;
import vdaoengine.utils.VSimpleProcedures;

public class ProcessTxtData {
	
	public static int numberOfDigitsToKeep = 3;

	public static void main(String[] args) {
		try{
			
			//CompileAandSTables("C:/Datas/BIODICA/work/OVCA_ICA/","OVCA_ica"); System.exit(0);
			//CompileAandSTablesAllResults("C:/Datas/MOSAIC/analysis/ica/metaanalysis/MOSAIC_ASP14/","Rescue_ica"); System.exit(0);
			//CompileAandSTablesAllResults("C:/Datas/MOSAIC/analysis/ica/metaanalysis/EMTAB/","EMTAB_ica"); System.exit(0);
			//CompileAandSTablesAllResults("C:/Datas/MOSAIC/analysis/ica/metaanalysis/ALLDATA_OLIVIER/","scmeta_ica"); System.exit(0);
			
			VDatReadWrite.writeNumberOfColumnsRows = false;
			VDatReadWrite.useQuotesEverywhere = false;

			
			if(args.length>0){
				String fn = args[0];
				VDataTable vt = null;
				if(!args[0].equals("-assembleICAresults")){
					vt = VDatReadWrite.LoadFromSimpleDatFile(fn, true, "\t", true);
					VSimpleProcedures.findAllNumericalColumns(vt);
				}
				for(int i=0;i<args.length;i++){
					if(args[i].equals("-logx1")){
						System.out.println("Convert to logarithmic scale x'<-log10(x+1)...");
						for(int k=0;k<vt.rowCount;k++)
							for(int j=0;j<vt.colCount;j++)if(vt.fieldTypes[j]==vt.NUMERICAL){
								float v = Float.parseFloat(vt.stringTable[k][j]);
								v = (float)Math.log10((double)v+1);
							      String fs = "#.";
							      for(int l=0;l<numberOfDigitsToKeep;l++)
							        fs+="#";
							      DecimalFormat df = new DecimalFormat(fs);
							      String s = df.format(v);
								vt.stringTable[k][j] = s;
							}
					}
					if(args[i].equals("-log2x1")){
						System.out.println("Convert to logarithmic scale x'<-log2(x+1)...");
						for(int k=0;k<vt.rowCount;k++)
							for(int j=0;j<vt.colCount;j++)if(vt.fieldTypes[j]==vt.NUMERICAL){
								float v = Float.parseFloat(vt.stringTable[k][j]);
								v = (float)(Math.log10((double)v+1)/Math.log10(2));
							      String fs = "#.";
							      for(int l=0;l<numberOfDigitsToKeep;l++)
							        fs+="#";
							      DecimalFormat df = new DecimalFormat(fs);
							      String s = df.format(v);
								vt.stringTable[k][j] = s;
							}
					}
					if(args[i].equals("-filtersum")){
						System.out.println("Filtering by sum of values in a row...");
						float val = Float.parseFloat(args[i+1]);
					}
					if(args[i].equals("-center")){
						System.out.println("Centering...");
						vt = VSimpleProcedures.normalizeVDat(vt, true, false);
					}
					if(args[i].equals("-doublecenter")){
						System.out.println("Double centering...");
						VDataSet vd = VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(vt, -1);
						float mas[][] = TableUtils.doubleCenterMatrix(vd.massif);
						for(int l=0;l<vt.rowCount;l++){
							int k=0;
							for(int s=0;s<vt.colCount;s++){
								if(vt.fieldTypes[s]==vt.NUMERICAL){
									vt.stringTable[l][s] = ""+mas[l][k];
									k++;
								}
							}
						}
					}
					if(args[i].equals("-selectrowsvar")){
						System.out.println("Selecting most variable "+args[i+1]+" rows...");
						vt = TableUtils.filterByVariation(vt, Integer.parseInt(args[i+1]), false);
					}
					if(args[i].equals("-selectuniquebyvar")){
						System.out.println("Selecting unique rows id based on variance...");
						vt = TableUtils.selectUniqueRowsIdsByVariance(vt, false);
						System.out.println(vt.rowCount+" were selected.");
					}
					if(args[i].equals("-selectrowsfromfile")){
						String filen = args[i+1];
						System.out.println("Selecting rows from file..."+filen);
						Vector<String> lines = Utils.loadStringListFromFile(filen);
						Vector<String> rows = new Vector<String>();
						rows.add(vt.fieldNames[0]);
						for(String s: lines){
							StringTokenizer st = new StringTokenizer(s,"\t");
							String row = st.nextToken();
							rows.add(row);
						}
						vt.makePrimaryHash(vt.fieldNames[0]);
						vt = VSimpleProcedures.selectRowsFromList(vt, rows);
					}
					if(args[i].equals("-ordercolumnspc")){
						int k = Integer.parseInt(args[i+1]);
						System.out.println("Ordering samples accordingly to PC1...");
						PCAMethod pca = new PCAMethod();
						VDataSet vd = VSimpleProcedures.SimplyPreparedDataset(vt, -1);
						pca.setDataSet(vd);
						pca.calcBasis(k);
						float contrs[] = new float[vd.coordCount];
						for(int j=0;j<vd.coordCount;j++){
							contrs[j] = (float)pca.getBasis().basis[k-1][j];
						}
						int inds[] = Algorithms.SortMass(contrs);
						Vector<String> numericalColumns = new Vector<String>();
						for(int j=0;j<vd.coordCount;j++){
							//System.out.println(vt.fieldNames[vd.selector.selectedColumns[inds[j]]]+"\t"+pca.getBasis().basis[k-1][inds[j]]);
							numericalColumns.add(vt.fieldNames[vd.selector.selectedColumns[inds[j]]]);
						}
						int l=0;
						Vector<String> columns = new Vector<String>();
						for(int j=0;j<vt.colCount;j++){
							if(vt.fieldTypes[j]==vt.STRING) columns.add(vt.fieldNames[j]);
							else{
								columns.add(numericalColumns.get(l)); l++;
							}
						}
						vt = VSimpleProcedures.SelectColumns(vt, columns);
					}
					if(args[i].equals("-selectcolumnsfromfile")){
						String filen = args[i+1];
						System.out.println("Selecting samples from file..."+filen);						
						Vector<String> lines = Utils.loadStringListFromFile(filen);
						Vector<String> cols = new Vector<String>();
						cols.add(vt.fieldNames[0]);
						for(String s: lines){
							StringTokenizer st = new StringTokenizer(s,"\t");
							String col = st.nextToken();
							if(vt.fieldNumByName(col)!=-1)
								cols.add(col);
						}
						vt = VSimpleProcedures.SelectColumns(vt, cols);
						System.out.println(vt.fieldNames[0]+"..."+vt.fieldNames[1]+"..."+vt.fieldNames[2]+"...");
					}
					if(args[i].equals("-decompose")){
						String sampleFile = args[i+1];
						StringTokenizer st = new StringTokenizer(sampleFile,"#");
						sampleFile = st.nextToken();
						String field = st.nextToken();
						System.out.println("Decomposing from... "+sampleFile);
						HashMap<String, Vector<String>> groupSample = new HashMap<String, Vector<String>>();
						Vector<String> allsamples = new Vector<String>();
						VDataTable samples = VDatReadWrite.LoadFromSimpleDatFile(sampleFile, true, "\t");
						for(int k=0;k<samples.rowCount;k++){
							String sample = samples.stringTable[k][0];
							String group = samples.stringTable[k][samples.fieldNumByName(field)];
							Vector<String> ss = groupSample.get(group);
							if(ss==null) ss = new Vector<String>();
							ss.add(sample);
							groupSample.put(group, ss);
							if(!allsamples.contains(sample)) 
								allsamples.add(sample);
						}
						Set<String> keys = groupSample.keySet();
						for(String group: keys){
							Vector<String> cols = groupSample.get(group);
							cols.insertElementAt(vt.fieldNames[0], 0);
							VDataTable vtg = VSimpleProcedures.SelectColumns(vt, cols);
							System.out.println("Saving "+group+".txt ...");
							VDatReadWrite.saveToSimpleDatFile(vtg, group+".txt");
						}
						System.out.println("Saving all groups ...");
						allsamples.insertElementAt(vt.fieldNames[0], 0);
						VDataTable vta = VSimpleProcedures.SelectColumns(vt, allsamples);
						VDatReadWrite.saveToSimpleDatFile(vta, fn.substring(0, fn.length()-4)+"_"+field+".txt");
					}
					if(args[i].equals("-collapsevar")){
						String sampleFile = args[i+1];
						StringTokenizer st = new StringTokenizer(sampleFile,"#");
						sampleFile = st.nextToken();
						String field = st.nextToken();
						int type = Integer.parseInt(st.nextToken());
						System.out.println("Collapsing variance from... "+sampleFile);
						HashMap<String, Vector<String>> groupSample = new HashMap<String, Vector<String>>();
						Vector<String> allsamples = new Vector<String>();
						VDataTable samples = VDatReadWrite.LoadFromSimpleDatFile(sampleFile, true, "\t");
						for(int k=0;k<samples.rowCount;k++){
							String sample = samples.stringTable[k][0];
							String group = samples.stringTable[k][samples.fieldNumByName(field)];
							Vector<String> ss = groupSample.get(group);
							if(ss==null) ss = new Vector<String>();
							ss.add(sample);
							groupSample.put(group, ss);
							if(!allsamples.contains(sample)) 
								allsamples.add(sample);
						}
						Set<String> keys = groupSample.keySet();
						Vector<String> keysv = new Vector<String>();
						for(String group: keys) keysv.add(group);
						Collections.sort(keysv);
						FileWriter fw = new FileWriter(fn.substring(0, fn.length()-4)+"_"+field+".txt");
						fw.write(vt.fieldNames[0]+"\t"); for(String group: keysv){ fw.write(field+"_"+group+"_var"+type+"\t"); } fw.write("\n");
						for(int row=0;row<vt.rowCount;row++){
							fw.write(vt.stringTable[row][0]+"\t");
						for(String group: keysv){
							Vector<String> cols = groupSample.get(group);
							//System.out.println(group+"\t"+vt.stringTable[row][0]);
							float f[] = new float[cols.size()];
							for(int kk=0;kk<f.length;kk++){
								f[kk] = Float.parseFloat(vt.stringTable[row][vt.fieldNumByName(cols.get(kk))]);
								//System.out.print(f[kk]+"\t");
							}
							//System.out.println();
							float val = 0f;
							if(type==0) val = VSimpleFunctions.calcStandardDeviation(f);
							if(type==2) val = VSimpleFunctions.calcStandardDeviationBiggerThan(f, 1e-6f);
							if(type==3) {
								float stdv = VSimpleFunctions.calcStandardDeviation(f);
								float mean = VSimpleFunctions.calcMean(f);
								val = stdv/Math.abs(mean+0.001f);
								val = val*val;
							}
							if(type==4) {
								val = VSimpleFunctions.calcMean(f);
							}
							fw.write(val+"\t");
							//System.out.println(val);
						}
						fw.write("\n");
						}
						fw.close();
						
					}
					if(args[i].equals("-mergefile")){
						String filen = args[i+1];
						System.out.println("Merging with file..."+filen);
						VDataTable vt1 = VDatReadWrite.LoadFromSimpleDatFile(filen, true, "\t");
						vt = VSimpleProcedures.MergeTables(vt, vt.fieldNames[0], vt1, vt1.fieldNames[0], "_");
					}
					if(args[i].equals("-transpose")){
						System.out.println("Transposing...");
						vt = vt.transposeTable(vt.fieldNames[0]);
					}
					if(args[i].equals("-prepare4ICA")){
						System.out.println("Preparing for ICA...");
						VDatReadWrite.writeNumberOfColumnsRows = false;
						VDatReadWrite.useQuotesEverywhere = false;
						VDatReadWrite.saveToSimpleDatFile(vt, fn.substring(0, fn.length()-4)+"_ica.txt");						
						prepareTable4ICA(fn.substring(0, fn.length()-4)+"_ica.txt");
					}
					if(args[i].equals("-savetxt")){
						System.out.println("Saving to..."+fn.substring(0, fn.length()-4)+"_copy.txt");
						VDatReadWrite.useQuotesEverywhere = false;
						VDatReadWrite.numberOfDigitsToKeep = 3;
						VDatReadWrite.writeNumberOfColumnsRows = false;
						VDatReadWrite.saveToSimpleDatFile(vt, fn.substring(0, fn.length()-4)+"_copy.txt");
					}
					if(args[i].equals("-assembleICAresults")){
						String fff = args[i+1];
						StringTokenizer st = new StringTokenizer(fff,";");
						String folder = st.nextToken();
						String prefix = st.nextToken();
						System.out.println("Assembling ICA results in "+folder+"...");
						CompileAandSTables(folder,prefix);
					}
				}
				if(vt!=null){
					System.out.println("Saving dat file to "+fn+".dat"+"...");
					VDatReadWrite.numberOfDigitsToKeep = 3;
					VDatReadWrite.saveToVDatFile(vt, fn+".dat");
				}
			}else{
				System.out.println("Please indicate the name of the tab-limited file with the first line containint column names to convert to ViDaExpert .dat format!");
				System.out.println("Options: [-center] [-doublecenter] [-selectrowsvar n] [-selectrowsfromfile filename] [-selectcolumnsfromfile filename] [-selectuniquebyvar] [-mergefile filename] [-prepare4ICA] [-ordercolumnspc1] [-savetxt] [-transpose] [-logx1] [-decompose sampleFileName#fieldName] [-collapsevar sampleFileName#fieldName#type]");
				System.out.println("Notes: for collapsevar type is 0 - standard variance, 1 - trimmed top 5%, 2 - exclude zeros, 3 - coeffvar^2, 4 - standard mean");
				
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}

	}
	
	public static void callFromMatlab(String args_string){
		System.out.println("Calling ProcessTxtData with arguments "+args_string);
		Vector<String> varg = new Vector<String>();
		StringTokenizer st = new StringTokenizer(args_string,"%");
		while(st.hasMoreTokens()){
			varg.add(st.nextToken());
		}
		String args[] = new String[varg.size()];
		for(int i=0;i<varg.size();i++) args[i] = varg.get(i);
		main(args);
	}
	
	  public static void prepareTable4ICA(String fname) throws Exception{
		  VDataTable vt = VDatReadWrite.LoadFromSimpleDatFile(fname, true, "\t");
		  String fn = fname.substring(0, fname.length()-4);
		  FileWriter fw = new FileWriter(fn+"_numerical.txt");
		  for(int i=0;i<vt.rowCount;i++){
			  for(int j=1;j<vt.colCount;j++)
				  fw.write(vt.stringTable[i][j]+"\t");
			  fw.write("\n");
		  }
		  fw.close();
		  FileWriter fws = new FileWriter(fn+"_samples.txt");
		  for(int j=1;j<vt.colCount;j++)
			  fws.write(vt.fieldNames[j]+"\t");
		  fws.close();
		  FileWriter fwids = new FileWriter(fn+"_ids.txt");
		  for(int i=0;i<vt.rowCount;i++){
			  fwids.write(vt.stringTable[i][0]+"\n");
		  }
		  fwids.close();
	  }
	  
	  public static void CompileAandSTables(String folder, String prefix) throws Exception{
		  CompileAandSTables(folder,prefix,-1);
	  }
	  
		public static void CompileAandSTables(String folder, String prefix, int numberOfComponentsToTake) throws Exception{

			Locale.setDefault(new Locale("en", "US"));
			System.out.println("Forcing using point as decimal separator... "+Locale.getDefault().getCountry());
			
			String ending = "_numerical.txt.num";
			// If file with exact prefix exists then we just proceed 
			if(new File(folder+"A_"+prefix).exists()){
				ending = "";
			}else{
			// First, let us determine the right file, with _numerical.txt_XX.num ending, we want to know XX
			File lf[] = new File(folder).listFiles();
			int maxnumcomp = 0;
			String selectedEnding = "";
			for(File f:lf){
				String fn = f.getName(); 
				if(fn.contains("_numerical.txt"))if(fn.endsWith(".num")) { 
					int k=fn.indexOf("_numerical.txt"); 
					ending = fn.substring(k, fn.length());
					//System.out.println(ending);
					int numcomp = Integer.parseInt(ending.substring(15,ending.length()-4));
					if(numcomp>maxnumcomp){
						maxnumcomp = numcomp;
						selectedEnding = ending;
					}
					if(numberOfComponentsToTake>0){
						if(numcomp==numberOfComponentsToTake)
							selectedEnding = ending;
					}
				}
			}
			
			/*if(!new File(folder+"A_"+prefix+ending).exists()){
				ending = "_ica"+ending;
			}*/
			ending = selectedEnding;
			}
						
			if(new File(folder+"A_"+prefix+ending).exists()){
			
			VDataTable vtA = VDatReadWrite.LoadFromSimpleDatFile(folder+"A_"+prefix+ending, false, "\t");
			System.out.println(vtA.rowCount+" samples in "+folder+"A_"+prefix+ending);
			VDataTable vtS = VDatReadWrite.LoadFromSimpleDatFile(folder+"S_"+prefix+ending, false, "\t");
			
			
			
			//Set<String> keys = sampleAnnotations.tableHashPrimary.keySet();
			//for(String sss: keys) System.out.println("Key '"+sss+"'");

			String prefixsid = prefix;
			if(ending.equals("")){
				int k = prefix.indexOf("_numerical.txt");
				prefixsid = prefix.substring(0, k);
			}
			
			String colString = Utils.loadString(folder+prefixsid+"_samples.txt");
			Vector<String> samples = new Vector<String>();
			StringTokenizer st = new StringTokenizer(colString,"\t");
			while(st.hasMoreTokens()){ 
				String ss = st.nextToken(); 
				if(!ss.trim().equals("")) samples.add(ss); 
			}
			System.out.println(samples.size()+" samples in "+folder+prefixsid+"_samples.txt");
			Vector<String> ids = Utils.loadStringListFromFile(folder+prefixsid+"_ids.txt");
			
			/*
			 * Normalize table of component projections, heavy tail always positive
			 */
			Vector<Boolean> swapped = new Vector<Boolean>();
			Vector<MetaGene> sMetagenes = new Vector<MetaGene>();
			for(int i=0;i<vtS.colCount;i++){
				String fn = vtS.fieldNames[i];
				Boolean swpd = false;
				System.out.println("Field "+fn);
				MetaGene mg = new MetaGene();
				for(int j=0;j<vtS.rowCount;j++){
					//if(j==100000*(int)(0.00001f*j))
					//	System.out.print(j+"\t");
					mg.add(ids.get(j), Float.parseFloat(vtS.stringTable[j][i]));
				}
				//System.out.println();
				//float positiveStd = mg.sidedStandardDeviation(+1);
				//float negativeStd = mg.sidedStandardDeviation(-1);
				float positiveSide = mg.sumOfWeightsAboveThreshold(+1,3f);
				float negativeSide = mg.sumOfWeightsAboveThreshold(-1,3f);
				if(positiveSide+negativeSide<1e-6){
					positiveSide = mg.sumOfWeightsAboveThreshold(+1,2f);
					negativeSide = mg.sumOfWeightsAboveThreshold(-1,2f);
					if(positiveSide+negativeSide<1e-6){
						positiveSide = mg.sumOfWeightsAboveThreshold(+1,1.5f);
						negativeSide = mg.sumOfWeightsAboveThreshold(-1,1.5f);
					}
				}
				
				System.out.println("positiveSide = "+positiveSide+", negativeSide = "+negativeSide);
				if(negativeSide>positiveSide){
					//System.out.println("TSPAN6 weight = "+mg.getWeight("TSPAN6"));
					mg.invertSignsOfWeights();
					//System.out.println("Flipped, TSPAN6 weight = "+mg.getWeight("TSPAN6"));
					swpd = true;
				}
				sMetagenes.add(mg);
				swapped.add(swpd);
			}

			System.out.println("Saving to "+folder+prefix+"_A.xls");
			
			
			FileWriter fwA = new FileWriter(folder+prefix+"_A.xls");
			fwA.write("SAMPLE\t"); for(int i=1;i<=vtA.colCount;i++) fwA.write("IC"+i+"\t"); fwA.write("\n");
			for(int i=0;i<vtA.rowCount;i++){
				String sample = samples.get(i);
				fwA.write(sample+"\t");
				for(int j=0;j<vtA.colCount;j++){
					float f = Float.parseFloat(vtA.stringTable[i][j]);
					if(swapped.get(j)) f = -f;
					DecimalFormat nf = new DecimalFormat("#.####");
					String vs = nf.format(f);
					fwA.write(vs+"\t");
				}
				fwA.write("\n");
			}
			fwA.close();
			
			System.out.println("Saving to "+folder+prefix+"_S.xls");
			FileWriter fwS = new FileWriter(folder+prefix+"_S.xls");
			fwS.write("PROBE\t"); for(int i=1;i<=vtS.colCount;i++) fwS.write("IC"+i+"\t"); fwS.write("\n");
			for(int i=0;i<vtS.rowCount;i++){
				String id = ids.get(i);
				if(i==100000*(int)(0.00001f*i))
					System.out.print(i+"\t");
				fwS.write(id+"\t");
				for(int j=0;j<vtS.colCount;j++){
					/*float f = Float.parseFloat(vtS.stringTable[i][j]);
					DecimalFormat nf = new DecimalFormat("#.###");
					String vs = nf.format(f);
					fwS.write(vs+"\t");
					if(idAnnotations!=null) fwSa.write(vs+"\t");*/
					float f = sMetagenes.get(j).getWeight(id);
					DecimalFormat nf = new DecimalFormat("#.####");
					String vs = nf.format(f);
					fwS.write(vs+"\t");
				}
				fwS.write("\n");
			}
			fwS.close();
			System.out.println();
			
			}else{
				System.out.println("Did not find "+folder+"A_"+prefix+ending);
			}
		}

		public static void CompileAandSTablesAllResults(String folder, String prefix) throws Exception{
			
			// First, let us determine the right file, with _numerical.txt_XX.num ending, we want to know XX
			File lf[] = new File(folder).listFiles();
			String ending = "_numerical.txt.num";
			Vector<String> processed = new Vector<String>();
			for(File f:lf){
				String fn = f.getName(); 
				if(fn.contains("_numerical.txt"))if(fn.endsWith(".num")) { 
					int k=fn.indexOf("_numerical.txt"); ending = fn.substring(k, fn.length());
					if(!processed.contains(prefix+ending)){
						System.out.println(prefix+ending);
						CompileAandSTables(folder,prefix+ending);
						processed.add(prefix+ending);
					}
				}
			}
			
		}

	  


}
