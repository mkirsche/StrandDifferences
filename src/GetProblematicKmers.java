import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.TreeSet;

public class GetProblematicKmers 
{
	
	static String tableFn = "";
	static String ofn = "";
	
	static int k = 6;
	static int occurrenceThreshold = 30;
	
	static int uniquePositionsThreshold = 3;
	
	static void usage()
	{
		System.out.println("Usage: java -cp src GetProblematicKmers [args]");
		System.out.println("  Example: java -cp src GetProblematicKmers table_file=table.txt out_file=kmers.txt");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  table_file   (String) - table of sites (with context) which cause strand bias issues");
		System.out.println("  out_file     (String) - file to record k-mers which seem to affect one strand but not the other");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  k                      (int)  [6] - length of k-mers to use");
		System.out.println("  occurrence_threshold   (int) [20] - number of times a k-mer must affect strand bias to be considered problematic");
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
				else if(key.equalsIgnoreCase("k")) { k = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("occurrence_threshold")) { occurrenceThreshold = Integer.parseInt(val); }
			}
		}
		
		if(tableFn.length() == 0 || ofn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		Scanner input = new Scanner(new FileInputStream(new File(tableFn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		
		String headerLine = input.nextLine();
		
		Table table = new Table(headerLine);
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			table.addRow(line);
		}
		
		TreeSet<KmerData> kmerData = new TreeSet<KmerData>();
		for(int i = 0; i<table.rows.size(); i++)
		{
			table.updateKmerData(i, kmerData);
		}
		
		out.println("KMER\tRC_KMER\tALT_KMER\tALT_RC_KMER\tCOUNT\tRC_COUNT\tSAMPLES\tPOSITIONS");
		for(KmerData kd : kmerData)
		{
			String s = kd.kmer;
			String alt = kd.altKmer;
			if(kd.count + kd.rcCount >= occurrenceThreshold || kd.uniquePositions.size() >= uniquePositionsThreshold)
			{
				String rc = reverseComplement(s);
				String rcAlt = reverseComplement(alt);
				String sampleList = "";
				String posList = "";
				HashSet<String> uniqueSamples = new HashSet<String>();
				HashSet<Integer> uniquePositions = new HashSet<Integer>();
				for(int i = 0; i<kd.samples.size(); i++)
				{
					if(!uniqueSamples.contains(kd.samples.get(i)))
					{
						if(sampleList.length() > 0)
						{
							sampleList += ", ";
						}
						sampleList += kd.samples.get(i);
					}
					if(!uniquePositions.contains(kd.positions.get(i)))
					{
						if(posList.length() > 0)
						{
							posList += ", ";
						}
						posList += kd.positions.get(i);
					}
					
					uniqueSamples.add(kd.samples.get(i));
					uniquePositions.add(kd.positions.get(i));
				}
				
				out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", s, rc, alt, rcAlt, kd.count, kd.rcCount, sampleList, posList);
			}
		}
		
		input.close();
		out.close();
	}
	
	static class KmerData implements Comparable<KmerData>
	{
		int count = 0;
		int rcCount = 0;
		String kmer;
		String altKmer;
		ArrayList<String> samples;
		ArrayList<Integer> positions;
		HashSet<Integer> uniquePositions;
		KmerData(String kmer, String altKmer, int pos, String sample)
		{
			count = 0;
			rcCount = 0;
			this.kmer = kmer;
			this.altKmer = altKmer;
			positions = new ArrayList<Integer>();
			samples = new ArrayList<String>();
			uniquePositions = new HashSet<Integer>();
			samples.add(sample);
			positions.add(pos);
			uniquePositions.add(pos);
		}
		
		void merge(KmerData kd)
		{
			count = count + kd.count;
			rcCount = rcCount + kd.rcCount;
			for(String sample : kd.samples)
			{
				samples.add(sample);
			}
			
			for(int pos : kd.positions)
			{
				positions.add(pos);
				uniquePositions.add(pos);
			}
		}
		
		@Override
		public int compareTo(KmerData o) {
			int res = kmer.compareTo(o.kmer);
			if(res != 0) 
			{
				return res;
			}
			return altKmer.compareTo(o.altKmer);
		}
	}
	
	static class Table
	{
		// Map of category names to which column they correspond to
		HashMap<String, Integer> categoryToIndex;
		
		// A list of data entries (rows in the table)
		ArrayList<String[]> rows;
		
		/*
		 * Parses the header line and initializes a table with those fields
		 */
		Table(String headerLine)
		{
			categoryToIndex = new HashMap<String, Integer>();
			String[] categories = headerLine.split("\t");
			for(int i = 0; i<categories.length; i++)
			{
				categoryToIndex.put(categories[i].toLowerCase(), i);
			}
			
			rows = new ArrayList<String[]>();
		}
		
		/*
		 * Adds a row to the table based on a TSV row
		 */
		void addRow(String line)
		{
			// Ignore empty lines
			if(line.length() == 0)
			{
				return;
			}
			String[] tokens = line.split("\t");
			
			rows.add(tokens);
		}
		
		/*
		 * Gets a particular field's value in a given row
		 */
		String getValue(int rowIndex, String category)
		{
			return rows.get(rowIndex)[categoryToIndex.get(category.toLowerCase())];
		}
		
		void updateKmerData(int rowIndex, TreeSet<KmerData> kmerData)
		{
			double plusMaf = Double.parseDouble(getValue(rowIndex, "PLUS_MAF"));
			double minusMaf = Double.parseDouble(getValue(rowIndex, "MINUS_MAF"));
			
			boolean plusStrand = plusMaf > minusMaf;
			String context = getValue(rowIndex, plusStrand ? "REF_CONTEXT" : "REF_CONTEXT_RC");
			
			char altChar = getAlt(rowIndex, plusStrand);
			
			for(int i = 0; i+k <= context.length(); i++)
			{
				String kmer = context.substring(i, i+k);
				
				int capIndex = -1;
				for(int j = 0; j<kmer.length(); j++)
				{
					if(Character.isUpperCase(kmer.charAt(j)))
					{
						capIndex = j;
						break;
					}
				}
				
				if(capIndex == -1)
				{
					// All lowercase so doesn't include the variant site
					continue;
				}
				
				kmer = kmer.toUpperCase();
				
				String altKmer = kmer.substring(0, capIndex) + Character.toUpperCase(altChar) + kmer.substring(1 + capIndex); 
				
				String rcKmer = reverseComplement(kmer);
				String altRcKmer = reverseComplement(altKmer);
				
				boolean usingMainKmer = kmer.compareTo(rcKmer) <= 0;
				String key = usingMainKmer ? kmer : rcKmer;
				String altKey = usingMainKmer ? altKmer : altRcKmer;
				KmerData cur = new KmerData(key, altKey, Integer.parseInt(getValue(rowIndex, "POS")), getValue(rowIndex, "SAMPLE"));
				if(!kmerData.contains(cur))
				{
					kmerData.add(cur);
				}
				else
				{
					kmerData.floor(cur).merge(cur);
				}
				
				KmerData toUpdate = kmerData.floor(cur);
				if(usingMainKmer)
				{
					toUpdate.count++;
				}
				else
				{
					toUpdate.rcCount++;
				}
			}
		}
		
		char getAlt(int rowIndex, boolean plusStrand)
		{
			int refVal =  GetStrandDifferences.charToInt(getValue(rowIndex, "REF").charAt(0));
			String[] alleleFreqs = getValue(rowIndex, plusStrand ? "PLUS_STRAND_FREQUENCIES" : "MINUS_STRAND_FREQUENCIES").split(",");
			int maxi = -1;
			int maxval = 0;
			for(int i = 0; i<4; i++)
			{
				int count = Integer.parseInt(alleleFreqs[i]);
				if(i == refVal)
				{
					continue;
				}
				if(maxi == -1 || count > maxval)
				{
					maxi = i;
					maxval = count;
				}
			}
			return GetStrandDifferences.intToChar(maxi);
		}
		
	}
	
	static String reverseComplement(String s)
	{
		int n = s.length();
		char[] res = new char[n];
		for(int i = 0; i<n; i++)
		{
			char c = s.charAt(n - 1 - i);
			res[i] = c;
			if(c == 'A' || c == 'a') res[i] = (char)(c + ('T' - 'A'));
			else if(c == 'C' || c == 'c') res[i] = (char)(c + ('G' - 'C'));
			else if(c == 'G' || c == 'g') res[i] = (char)(c + ('C' - 'G'));
			else if(c == 'T' || c == 't') res[i] = (char)(c + ('A' - 'T'));
		}
		return new String(res);
	}

}
