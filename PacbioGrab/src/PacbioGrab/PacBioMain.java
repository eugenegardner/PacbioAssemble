package PacbioGrab;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.ggf.drmaa.DrmaaException;

import MELT.utilities.FastqWriter;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class PacBioMain {

	public static void main(String[] args) throws IOException, InterruptedException, DrmaaException {

		if (args.length != 6 || args[0].equals("-h") || args[0].equals("--help") || args[0].equals("-help")) {
			
			System.out.println("Pacbio.jar usage:");
			System.out.println("java -Xmx2g -jar Pacbio.jar <reference.fa> <human_aligned.sorted.bam> <L1_aligned.sorted.bam> <raw_movies.bam> <working_dir/> <fofn.txt>");
			System.exit(1);
			
		}
		
		//Check to make sure java is 1.8 and blasr is executable
		if (checkExes() == false) {
			System.exit(1);
		}
		
		File fastaIndex = new File(args[0] + ".fai");
		File pacBAM = new File(args[1]);
		File pacBAI = new File(args[1] + ".bai");
		File lineBAM = new File(args[2]);
		File lineBAI = new File(args[2] + ".bai");
		File rawBAM = new File(args[3]);
		String workingDir = args[4];
		setDrmaaLoc(); //Make sure Drmaa is included in path.
		
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		
		SamReader pacReader = readerFactory.open(SamInputResource.of(pacBAM).index(pacBAI));	
		BufferedReader indexReader = new BufferedReader(new FileReader(fastaIndex));
		SamReader lineReader = readerFactory.open(SamInputResource.of(lineBAM).index(lineBAI));
		SamReader rawReader = readerFactory.open(SamInputResource.of(rawBAM));
		
		String line;
		String data[];
		
		File wd = new File(workingDir);
		if (!wd.exists()) {
			wd.mkdir();
		}
		new File(workingDir + "/LINE1_alignments/").mkdir();
		new File(workingDir + "/FASTQ_Files/").mkdir();
		new File(workingDir + "/blasr_alignments/").mkdir();
		new File(workingDir + "/White_Files/").mkdir();
		new File(workingDir + "/tmp/").mkdir();
		
		
		Map<String, HashMap<String, Integer>> numTot = new HashMap<String, HashMap<String, Integer>>();
		Map<String, FastqWriter> fastqWriters = new HashMap<String, FastqWriter>();
		Map<String, SAMFileWriter> lineWriters = new HashMap<String, SAMFileWriter>();
		Map<String, SAMFileWriter> rawWriters = new HashMap<String, SAMFileWriter>();
		Map<String, BufferedWriter> whiteWriters = new HashMap<String, BufferedWriter>();
		
		SAMFileWriterFactory samFactory = new SAMFileWriterFactory().setCreateIndex(true).setTempDirectory(new File(workingDir + "/tmp/"));
		
		while ((line = indexReader.readLine()) != null) {
			
			data = line.split("\t");
			String name = data[0];
			new File(workingDir + "/" + name).mkdir();
			
			// Grab reads greater than 1500
			
			HashMap<String, Integer> readNames = PacBioMethods.getReads(pacReader, name);
			fastqWriters.put(name, new FastqWriter(new BufferedWriter(new FileWriter(new File(workingDir + "/FASTQ_Files/" + name + ".fq")))));
			lineWriters.put(name, samFactory.makeBAMWriter(lineReader.getFileHeader(), true, new File(workingDir + "/LINE1_alignments/" + name + ".sorted.bam")));
			rawWriters.put(name, samFactory.makeBAMWriter(rawReader.getFileHeader(), false, new File(workingDir + "/LINE1_alignments/" + name + ".raw.bam")));
			whiteWriters.put(name, new BufferedWriter(new FileWriter(new File(workingDir + "/White_Files/" + name + ".txt"))));
			numTot.put(name, readNames);
			
		}

		indexReader.close();
		pacReader.close();
		//Now go through the L1 bam and put reads in proper categories:
		PacBioMethods.printLINE(lineReader, rawReader, numTot, fastqWriters, lineWriters, whiteWriters, rawWriters, workingDir);
		//Close all file handles:
		for (Map.Entry<String, HashMap<String, Integer>> entry : numTot.entrySet()) {
			fastqWriters.get(entry.getKey()).closeBuffer();
			lineWriters.get(entry.getKey()).close();
			whiteWriters.get(entry.getKey()).close();
			rawWriters.get(entry.getKey()).close();
		}
		lineReader.close();
		System.exit(0);
		PacBioMethods.runConsensus(workingDir, numTot, new File(args[4]));

		
		// Now analyze the data for each ConsenseTools run
		
		System.out.println("CHR\tPOS\tNUM READS\tNUM SEQ");
		
		for (Map.Entry<String, HashMap<String, Integer>> entry : numTot.entrySet()) {
			
			String name = entry.getKey();
			int totalReads = entry.getValue().size();
			
			Pattern namePatt = Pattern.compile("([0-9X]{1,2})_(\\d+)");
			Matcher nameMatch = namePatt.matcher(name);
			
			String chr = null;
			String pos = null;
			if (nameMatch.matches()) {
				chr = nameMatch.group(1);
				pos = nameMatch.group(2);
			}
			if (totalReads <= 10) {
				
				System.out.println(chr + "\t" + pos + "\t" + totalReads + "\t0");
				
			} else {
				
				Map<Integer, Integer> seqs = PacBioMethods.analyzeConsenseFasta(workingDir, name);
				
				int numSeqs = seqs.size();
				
				if (numSeqs == 0) {
					
					System.out.println(chr + "\t" + pos + "\t" + totalReads + "\t" + numSeqs);
					
				} else {
				
					System.out.print(chr + "\t" + pos + "\t" + totalReads + "\t" + numSeqs + "\t");
				
					for (Map.Entry<Integer, Integer> total : seqs.entrySet()) {
					
						if (total.getKey() == numSeqs) {
						
							System.out.println(total.getKey() +"\t" + total.getValue());
						
						} else {
						
							System.out.print(total.getKey() + "\t" + total.getValue() + "\t");
						
						}
					
					}
				
				}
			
			}
			
		}
		
	}
	
	private static void setDrmaaLoc () {
		
		try {
			Field usrPathsField = ClassLoader.class.getDeclaredField("usr_paths");
			usrPathsField.setAccessible(true);
			final String[] paths = (String[])usrPathsField.get(null);
			final String[] newPaths = Arrays.copyOf(paths, paths.length+1);
			newPaths[newPaths.length-1] = "/usr/local/packages/sge-root/lib/lx24-amd64/";
			usrPathsField.set(null, newPaths);
		} catch (SecurityException e) {
			e.printStackTrace();
		} catch (NoSuchFieldException e) {
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		
	}
	private static boolean checkExes() throws IOException, InterruptedException {
		
		boolean javaWorks = false;
		boolean blasrWorks = false;
		Runtime rt = Runtime.getRuntime();
		Process pr = rt.exec("java -version");
		pr.waitFor();
		BufferedReader javaChecker = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
		String line;
		while ((line = javaChecker.readLine()) != null) {

			if (line.matches("java version \"1.8\\S+")) {
				javaWorks = true;
				break;
			}
			
		}
		javaChecker.close();
		String path = System.getenv("PATH");
		String paths[] = path.split(":");
		for (String p : paths) {
			File blasr = new File(p + "/blasr");
			if (blasr.exists()) {
				blasrWorks = true;
			}
		}
		
		boolean works = true;
		if (javaWorks == false) {
			System.out.println("java does not appear to be the correct version, needs java 1.8 (try 'use java1.8')");
			works = false;
		}
		if (blasrWorks == false) {
			System.out.println("blasr does not appear to be in PATH (try 'use blasr')");
			works = false;
		}
		return works;
	}

}
