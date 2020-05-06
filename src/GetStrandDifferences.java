import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;

public class GetStrandDifferences
{
	static String mpileupFn = "", ofn = "", genomeFn = "";
	
	static int maxLen = 31000;
	static int minDepth = 30;
	
	static int contextLength = 10;
	
	static double mafRatio = 2.0;
	static double minMaf = 0.15;
	
	// These are used when incorporating gene annotations
	static HashMap<String, String> genome;
	
	static void usage()
	{
		System.out.println("Usage: java -cp src GetStrandDifferences [args]");
		System.out.println("  Example: java -cp src GetStrandDifferences mpileup_file=mpileup.txt out_file=differences.txt genome_file=genome.fa");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  mpileup_file (String) - mpileup file");
		System.out.println("  out_file     (String) - file to record positions with strand differences");
		System.out.println("  genome_file  (String) - path to genome");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  min_depth (int)   [30]   - the minimum unambiguous depth that must be present on each strand for a position to be highlighted");
		System.out.println("  context   (int)   [10]   - the number of bases to report on either side of highlighted sites");
		System.out.println("  maf_ratio (float) [2.0]  - the minimum ratio of MAFs across strands needed to highlight a site");
		System.out.println("  min_maf   (float) [0.15] - the minimum MAF on the more frequent strand needed to highlight a site");
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
				if(key.equalsIgnoreCase("mpileup_file")) { mpileupFn = val; }
				else if(key.equalsIgnoreCase("out_file")) { ofn = val; } 
				else if(key.equalsIgnoreCase("genome_file")) { genomeFn = val; } 
				else if(key.equals("min_depth")) { minDepth = Integer.parseInt(val); }
				else if(key.equals("context")) { contextLength = Integer.parseInt(val); }
				else if(key.equals("maf_ratio")) { mafRatio = Double.parseDouble(val); }
				else if(key.equals("min_maf")) { minMaf = Double.parseDouble(val); }
			}
		}
		
		if(mpileupFn.length() == 0 || ofn.length() == 0)
		{
			usage();
			System.exit(1);
		}
		
		if(genomeFn.length() == 0)
		{
			usage();
			System.exit(1);;
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		genome = new HashMap<String, String>();

		// Read in genome
		Scanner genomeInput = new Scanner(new FileInputStream(new File(genomeFn)));

		StringBuilder seq = new StringBuilder("");
		String refName = "";

		while(genomeInput.hasNext())
		{
			String line = genomeInput.nextLine();
			if(line.startsWith(">"))
			{
				String name = line.split(" ")[0].substring(1);
				if(refName.length() > 0)
				{
					// add last contig
					genome.put(refName, seq.toString());
					seq = new StringBuilder("");
				}
				refName = name;
			}
			else
			{
				seq.append(line);
			}
		}
		if(refName.length() > 0)
		{
			// add last contig
			genome.put(refName, seq.toString());
		}
		genomeInput.close();
		
		findDifferences(mpileupFn, ofn);
		
	}
	
	/*
	 * Processes an mpileup file and highlight sites with major strand differences
	 */
	static void findDifferences(String mpileupFn, String ofn) throws Exception
	{
		Mpileup mp = new Mpileup(mpileupFn);
		PrintWriter out = new PrintWriter(new File(ofn));
		out.println("CHR\tPOS\tREF\tPLUS_STRAND_FREQUENCIES\tMINUS_STRAND_FREQUENCIES\tPLUS_MAF\tMINUS_MAF\tREF_CONTEXT\tREF_CONTEXT_RC");
		for(String chrName : mp.allFrequencies.keySet())
		{
			int[][][] counts = mp.allFrequencies.get(chrName);
			for(int i = 0; i<counts.length; i++)
			{
				int[] plusCounts = counts[i][1];
				int[] minusCounts = counts[i][2];
				
				int unambigPlusCov = plusCounts[0] + plusCounts[1] + plusCounts[2] + plusCounts[3];
				int unambigMinusCov = minusCounts[0] + minusCounts[1] + minusCounts[2] + minusCounts[3];
				
				if(unambigPlusCov < minDepth || unambigMinusCov < minDepth)
				{
					continue;
				}
				
				int maxPlus = -1, maxMinus = -1;
				
				char refChar = genome.get(chrName).charAt(i);
				int refVal = charToInt(refChar);
				
				for(int j = 0; j<4; j++)
				{
					if(j == refVal) continue;
					if(maxPlus == -1 || plusCounts[j] > plusCounts[maxPlus]) maxPlus = j;
					if(maxMinus == -1 || minusCounts[j] > minusCounts[maxMinus]) maxMinus = j;
				}
				
				double plusMaf = 1.0 * plusCounts[maxPlus] / unambigPlusCov;
				double minusMaf = 1.0 * minusCounts[maxMinus] / unambigMinusCov;
				
				double higherMaf = Math.max(plusMaf, minusMaf), lowerMaf = Math.min(plusMaf, minusMaf);
								
				int contextStart = Math.max(0, i - contextLength);
				int contextEnd = Math.min(i + contextLength + 1, genome.get(chrName).length());
				
				char[] context = genome.get(chrName).substring(contextStart, contextEnd).toLowerCase().toCharArray();
				context[i - contextStart] += 'A' - 'a';
				
				char[] revComp = new char[context.length];
				for(int j = 0; j<context.length; j++)
				{
					char c = context[context.length - 1 - j];
					if(c == 'A') revComp[j] = 'T';
					else if(c == 'C') revComp[j] = 'G';
					else if(c == 'G') revComp[j] = 'C';
					else if(c == 'T') revComp[j] = 'A';
					else if(c == 'a') revComp[j] = 't';
					else if(c == 'c') revComp[j] = 'g';
					else if(c == 'g') revComp[j] = 'c';
					else if(c == 't') revComp[j] = 'a';
					else revComp[j] = c;
				}
				
				if(higherMaf >= minMaf - 1e-9 && higherMaf >= lowerMaf * mafRatio - 1e-9)
				{
					System.out.println(higherMaf+" "+lowerMaf);
					out.printf("%s\t%s\t%s\t%d,%d,%d,%d,%d\t%d,%d,%d,%d,%d\t%.3f\t%.3f\t%s\t%s\n", chrName, i+1, refChar, 
							plusCounts[0], plusCounts[1], plusCounts[2], plusCounts[3], plusCounts[4],
							minusCounts[0], minusCounts[1], minusCounts[2], minusCounts[3], minusCounts[4],
							plusMaf, minusMaf, new String(context), new String(revComp));
				}
			}
		}
		out.close();
	}
	
	static class Mpileup
	{
		// Map chromosome name to an array of frequencies indexed by (position, strand, base)
		HashMap<String, int[][][]> allFrequencies;
		
		// The reference characters
		HashMap<String, char[]> genome;
		
		/*
		 * Take in an mpileup file and store the allele frequencies at each position
		 */
		Mpileup(String fn) throws Exception
		{
			Scanner input = new Scanner(new FileInputStream(new File(fn)));
			allFrequencies = new HashMap<String, int[][][]>();
			while(input.hasNext())
			{
				String line = input.nextLine();
				if(line.length() == 0 || line.startsWith("@"))
				{
					continue;
				}
				
				String[] tokens = line.split("\t");
				
				// Get chromosome, position, and ref allele
				String chrName = tokens[0];
				int refPos = Integer.parseInt(tokens[1]) - 1;
				char refChar = tokens[2].charAt(0);
				
				if(!allFrequencies.containsKey(chrName))
				{
					allFrequencies.put(chrName, new int[maxLen][3][6]);
				}
				
				// Fill the frequency array at this position
				int[][][] covArray = allFrequencies.get(chrName);
				covArray[refPos] = getAlleleFreqs(refChar, tokens[4]);
			}
			input.close();
		}
	}
	
	/*
	 * Gets the number of A/C/G/T/N's covering a position from an mpileup string
	 */
	static int[][] getAlleleFreqs(char refChar, String pileup)
	{
		int[][] res = new int[3][6];
		for(int i = 0; i<res.length; i++)
		{
			res[i] = new int[6];
		}
		for(int i = 0; i<pileup.length(); i++)
		{
			char c = pileup.charAt(i);
			
			// Exact match so use ref character
			if(c == '.' || c == ',')
			{
				res[0][charToInt(refChar)]++;
				if(c == '.')
				{
					res[1][charToInt(refChar)]++;
				}
				else
				{
					res[2][charToInt(refChar)]++;
				}
			}
			
			// Insertion or deletion after this base so ignore
			else if(c == '+' || c == '-')
			{
				int end = i;
				int length = 0;
				while(end+1 < pileup.length() && pileup.charAt(end+1) >= '0' && pileup.charAt(end+1)<= '9')
				{
					end++;
					length = length * 10 + pileup.charAt(end) - '0';
				}
				boolean capital = (pileup.charAt(end+1) >= 'A' && pileup.charAt(end+1) <= 'Z') || pileup.charAt(end+1) == '*';
				i = end + length;
				res[0][5]++;
				if(capital)
				{
					res[1][5]++;
				}
				else
				{
					res[2][5]++;
				}
			}
			
			else if(c == '*')
			{
				res[0][5]++;
				res[1][5]++;
			}
			
			else if(c == '#')
			{
				res[0][5]++;
				res[2][5]++;
			}
			
			// Last character indicator - ignore
			else if(c == '$')
			{
				continue;
			}
			
			// First character indicator - ignore
			else if(c == '^')
			{
				i++;
				continue;
			}
			
			// Mismatch or N so count this character after converting it to an integer
			else
			{
				int val = charToInt(c);
				if(val != -1)
				{
					res[0][charToInt(c)]++;
					
					if(Character.isUpperCase(c) || c == '>')
					{
						res[1][charToInt(c)]++;
					}
					if(Character.isLowerCase(c) || c == '<')
					{
						res[2][charToInt(c)]++;
					}
				}
			}
		}
		return res;
	}
	
	/*
	 * Converts a basepair charater to an integer index
	 */
	static int charToInt(char c)
	{
		if(c == 'a' || c == 'A') return 0;
		else if(c == 'c' || c == 'C') return 1;
		else if(c == 'g' || c == 'G') return 2;
		else if(c == 't' || c == 'T') return 3;
		else if(c == '>' || c == '<' || c == 'n' || c == 'N') return 4;
		else return -1;
	}

	static char intToChar(int val)
	{
		if(val == 0) return 'A';
		else if(val == 1) return 'C';
		else if(val == 2) return 'G';
		else if(val == 3) return 'T';
		else return 'N';
	}
}
