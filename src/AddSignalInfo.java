import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;

public class AddSignalInfo {
	static String tableFn = "", modelFn = "", ofn = "";
	static void usage()
	{
		System.out.println("Usage: java -cp src AddSignalInfo [args]");
		System.out.println("  Example: java -cp src AddSignalInfo table_file=kmers.txt model_file=model.txt out_file=kmers.signal.txt");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  table_file   (String) - table of k-mers which frequently occur in strand-specific variants");
		System.out.println("  out_file     (String) - file to write the updated table to");
		System.out.println("  model_file   (String) - file with mean and standard deviations of signal for each k-mer");
		System.out.println();
	}
	
	static void parseArgs(String[] args)
	{
		for(String s : args)
		{
			int equalsIdx = s.indexOf('=');
			if(equalsIdx == -1)
			{
				
			}
			else
			{
				String key = s.substring(0, equalsIdx);
				String val = s.substring(1 + equalsIdx);
				if(key.equalsIgnoreCase("table_file")) { tableFn = val; }
				else if(key.equalsIgnoreCase("out_file")) { ofn = val; } 
				else if(key.equalsIgnoreCase("model_file")) { modelFn = val; }
			}
		}
		
		if(tableFn.length() == 0 || ofn.length() == 0 || modelFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
public static void main(String[] args) throws Exception
{
	parseArgs(args);
	Scanner modelInput = new Scanner(new FileInputStream(new File(modelFn)));
	GetProblematicKmers.Table model = new GetProblematicKmers.Table(modelInput.nextLine());
	while(modelInput.hasNext())
	{
		model.addRow(modelInput.nextLine());
	}
	modelInput.close();
	
	HashMap<String, Double> kmerToMean = new HashMap<String, Double>();
	HashMap<String, Double> kmerToStdev = new HashMap<String, Double>();
	
	for(int i = 0; i<model.rows.size(); i++)
	{
		String kmer = model.getValue(i, "kmer");
		double mean = Double.parseDouble(model.getValue(i, "level_mean"));
		double stdev = Double.parseDouble(model.getValue(i, "level_stdv"));
		kmerToMean.put(kmer, mean);
		kmerToStdev.put(kmer, stdev);
	}
	
	Scanner input = new Scanner(new FileInputStream(new File(tableFn)));
	PrintWriter out = new PrintWriter(new File(ofn));
		
	String header = input.nextLine();
	
	GetProblematicKmers.Table table = new GetProblematicKmers.Table(header);

	out.println(header + 
			"\t" + "level_mean" + "\t" + "level_stdv" +
			"\t" + "alt_level_mean" + "\t" + "alt_level_stdv" + 
			"\t" + "rc_level_mean" + "\t" + "rc_level_stdv" +
			"\t" + "alt_rc_level_mean" + "\t" + "alt_rc_level_stdv");
	while(input.hasNext())
	{
		String line = input.nextLine();
		table.addRow(line);
		int rowIndex = table.rows.size() - 1;
		
		String forwardKmer = table.getValue(rowIndex, "kmer");
		String revKmer = table.getValue(rowIndex, "rc_kmer");
		String altKmer = table.getValue(rowIndex, "alt_kmer");
		String altRevKmer = table.getValue(rowIndex, "alt_rc_kmer");

		out.printf("%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", line,
				kmerToMean.get(forwardKmer), kmerToStdev.get(forwardKmer),
				kmerToMean.get(altKmer), kmerToStdev.get(altKmer), 
				kmerToMean.get(revKmer), kmerToStdev.get(revKmer),
				kmerToMean.get(altRevKmer), kmerToStdev.get(altRevKmer));
	}
	input.close();
	out.close();
}
}
