package PacbioGrab;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;

import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.JobTemplate;
import org.ggf.drmaa.Session;
import org.ggf.drmaa.SessionFactory;

import MELT.utilities.FastqWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class PacBioMethods {
	
	public PacBioMethods () {
		
		super();
		
	}
	
	public static HashMap<String, Integer> getReads(SamReader pacReader, String chr) {
		
		SAMRecordIterator itr = pacReader.queryOverlapping(chr, 450, 550);
		
		HashMap<String, Integer> readNames = new HashMap<String, Integer>();
		
		while (itr.hasNext()) {
			
			SAMRecord samRecord = itr.next();
			
			int readLength = samRecord.getReadLength();
			int mapQ = samRecord.getMappingQuality();
			String readName = samRecord.getReadName();
			
			if (readLength >= 1500 && mapQ > 0) {
				
				readNames.put(readName, readLength);
				
			}
		}
		
		itr.close();
		
		return readNames;
		
	}
	public static void printLINE(SamReader lineReader, SamReader rawReader, Map<String, HashMap<String, Integer>> numTot, Map<String,FastqWriter> fastqWriters, Map<String,SAMFileWriter> lineWriters, Map<String, BufferedWriter> whiteWriters, Map<String,SAMFileWriter> rawWriters, String workingDir) throws IOException {
		
		SAMRecordIterator lineItr = lineReader.iterator();
//		Pattern nPattRight = Pattern.compile("\\S*[ATCG](N+)");
//		Pattern nPattLeft = Pattern.compile("(N+)[ATCG]\\S*");
		while (lineItr.hasNext()) {
			
			SAMRecord samRecord = lineItr.next();
			if (samRecord.getReadUnmappedFlag()) {
				break;
			}
			for (Map.Entry<String, HashMap<String, Integer>> entry : numTot.entrySet()) {
				Map<String, Integer> readsToGrab = entry.getValue();
				if (readsToGrab.containsKey(samRecord.getReadName())) {
					
//					String readString = samRecord.getReadString();
//					String qualityString = samRecord.getBaseQualityString();
//					Matcher nLeft = nPattLeft.matcher(readString);
//					if (nLeft.matches()) {
//						readString = readString.substring(nLeft.group(1).length(), readString.length());
//						qualityString = qualityString.substring(nLeft.group(1).length(), qualityString.length());
//					}
//					Matcher nRight = nPattRight.matcher(readString);
//					if (nRight.matches()) {
//						readString = readString.substring(0, (readString.length() - nRight.group(1).length()));
//						qualityString = qualityString.substring(0, (qualityString.length() - nRight.group(1).length()));
//					}
					
//					fastqWriters.get(entry.getKey()).writeSequence(samRecord.getReadName(), readString, qualityString);
					lineWriters.get(entry.getKey()).addAlignment(samRecord);
					whiteWriters.get(entry.getKey()).write(samRecord.getReadName() + "\n");
					whiteWriters.get(entry.getKey()).flush();
					
				}
			
			}
				
		}
		
		lineItr.close();
		
		SAMRecordIterator rawItr = rawReader.iterator();
		while (rawItr.hasNext()) {
			
			SAMRecord samRecord = rawItr.next();
			for (Map.Entry<String, HashMap<String, Integer>> entry : numTot.entrySet()) {
				Map<String, Integer> readsToGrab = entry.getValue();
				if (readsToGrab.containsKey(samRecord.getReadName())) {
					
					rawWriters.get(entry.getKey()).addAlignment(samRecord);
					
				}
			}			
		}
		
	}
	public static void pbcr (String workingDir, String specFile, Map<String, LinePackage> numTot) throws IOException, InterruptedException, DrmaaException {
		
		List<String> jobIDs = new CopyOnWriteArrayList<String>();
		Session session = PacBioMethods.buildSession();
		
		for (Map.Entry<String, LinePackage> entry : numTot.entrySet()) {
			if (entry.getValue().getNumReads() > 10) {
				String name = entry.getKey();
				JobTemplate jt = buildJobTemplate(session, workingDir, name);
				String exec = "/home/eugene.gardner/EugeneTools/wgs-8.3rc1/Linux-amd64/bin/PBcR -l chr" + name + " -fastq " + entry.getValue().getFqFile().getAbsolutePath() + " -length 3000 -maxCoverage 3000 -t 16 -s " + specFile + " -partitions 64 -genomeSize 10000";
				jt.setRemoteCommand(exec);
				jt.setJobName("CONS_" + name);
				String id = session.runJob(jt);
				jobIDs.add(id);
			}
		}
		
		while (jobIDs.isEmpty() == false) {

			for (String jobID : jobIDs) {
				
				try {
					if (session.getJobProgramStatus(jobID) == Session.DONE) {
					
						jobIDs.remove(jobID);
						
					}
						
				} catch (DrmaaException e) {
					
				}
				
			}
			
		}
		
	}
	public static void runConsensus(String workingDir, Map<String, HashMap<String, Integer>> numTot, File fofn) throws DrmaaException {
		List<String> jobIDs = new CopyOnWriteArrayList<String>();
		Session session = PacBioMethods.buildSession();
		
		for (Map.Entry<String, HashMap<String, Integer>> entry : numTot.entrySet()) {
			if (entry.getValue().size() > 10) {
				String name = entry.getKey();
				JobTemplate jt = buildJobTemplate(session, workingDir, name);
				String exec = "/local/projects-t2/LIPAC/SMRTAnalysis/install/smrtanalysis_2.3.0.140936/analysis/bin/ConsensusTools.sh AmpliconAnalysis --noPhasing --noClustering -f " + fofn.getAbsolutePath() + " --whiteList=" + workingDir + "/White_Files/" + name + ".txt -r 500 --minPredictedAccuracy=0.99 -n 16";
				jt.setRemoteCommand(exec);
				jt.setJobName("CONS_" + name);
				String id = session.runJob(jt);
				jobIDs.add(id);
			}
		}
		
		while (jobIDs.isEmpty() == false) {

			for (String jobID : jobIDs) {
				
				try {
					if (session.getJobProgramStatus(jobID) == Session.DONE) {
					
						jobIDs.remove(jobID);
						
					}
						
				} catch (DrmaaException e) {
					
				}
				
			}
			
		}
		
	}
	public static Map<Integer, Integer> analyzePBCRFasta (String workingDir, String name) throws IOException, InterruptedException {
		
		File asmFile = new File(workingDir + "/" + name + "/chr" + name + "/9-terminator/asm.ctg.fasta");
		Map<Integer,Integer> seqs = new LinkedHashMap<Integer,Integer>();
		
		if (asmFile.exists()) {
			
			BufferedReader fastaReader = new BufferedReader(new FileReader(asmFile));
			//Run blasr on the fasta file:
			runBlasr(workingDir, name, asmFile.getAbsolutePath());
			
			int length = 0;
			int numSeqs = 0;
			boolean first = true;
			String line;
			
			while((line = fastaReader.readLine()) != null) {
			
				if (first == true && line.matches("^\\>\\S+")) {
	
					first = false;
					numSeqs++;
				
				} else if (first == false && line.matches("^\\>\\S+")) {
	
					seqs.put(numSeqs, length);
					length = 0;
					numSeqs++;
				
				} else {
				
					length += line.length();
				
				}
			
			}
	
			if (numSeqs > 0) { 
			
				seqs.put(numSeqs, length);
				fastaReader.close();
			
			}
			
		}
			
		return seqs;
		
		
	}
	public static Map<Integer, Integer> analyzeConsenseFasta (String workingDir, String name) throws IOException, InterruptedException {
		
		File asmFile = new File(workingDir + "/" + name + "/amplicon_analysis.fasta");
		Map<Integer,Integer> seqs = new LinkedHashMap<Integer,Integer>();
		
		if (asmFile.exists()) {
			
			BufferedReader fastaReader = new BufferedReader(new FileReader(asmFile));
			//Run blasr on the fasta file:
			runBlasr(workingDir, name, asmFile.getAbsolutePath());
			
			int length = 0;
			int numSeqs = 0;
			boolean first = true;
			String line;
			
			while((line = fastaReader.readLine()) != null) {
			
				if (first == true && line.matches("^\\>\\S+")) {
	
					first = false;
					numSeqs++;
				
				} else if (first == false && line.matches("^\\>\\S+")) {
	
					seqs.put(numSeqs, length);
					length = 0;
					numSeqs++;
				
				} else {
				
					length += line.length();
				
				}
			
			}
	
			if (numSeqs > 0) { 
			
				seqs.put(numSeqs, length);
				fastaReader.close();
			
			}
			
		}
			
		return seqs;
		
		
	}
	public static JobTemplate buildJobTemplate (Session session, String workingDir, String name) throws DrmaaException {
		
		JobTemplate jt = session.createJobTemplate();
		jt.setNativeSpecification("-q threaded.q -P sdevine-lab -l mem_free=8G -pe thread 16 -b y -shell yes -V -wd " + workingDir + "/" + name + "/");
		String gridError = workingDir + "/error/";
		String gridOutput = workingDir + "/output/";
		File gridErrorFile = new File(gridError);
		File gridOutputFile = new File(gridOutput);
		if (gridErrorFile.exists() == false) {
			gridErrorFile.mkdir();
		}
		if (gridOutputFile.exists() == false) {
			gridOutputFile.mkdir();
		}
		jt.setOutputPath(":" + gridOutput);
		jt.setErrorPath(":" + gridError);

		return jt;
		
	}
	public static Session buildSession() throws DrmaaException {
		
		Session session = SessionFactory.getFactory().getSession();
		session.init("");
		return session;
		
	}
	private static void runBlasr(String workingDir, String name, String faLoc) throws IOException, InterruptedException {
		
		Runtime rt = Runtime.getRuntime();
		Process pr = rt.exec("blasr -sa /local/aberdeen2rw/DEVINE/LAB_FOLDERS/EugeneAnalysis/me_refs/line_ref/LINE1.fa.sa -clipping soft -sam -out " + workingDir + "/blasr_alignments/" + name + ".sam " + faLoc + " /local/aberdeen2rw/DEVINE/LAB_FOLDERS/EugeneAnalysis/me_refs/line_ref/LINE1.fa");
		pr.waitFor();
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader blasrReader = readerFactory.open(SamInputResource.of(new File(workingDir + "/blasr_alignments/" + name + ".sam")));	
		SAMFileHeader header = blasrReader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter blasrWriter = new SAMFileWriterFactory().setCreateIndex(true).setTempDirectory(new File(workingDir + "/tmp/")).makeBAMWriter(header, false, new File(workingDir + "/blasr_alignments/" + name + ".sorted.bam"));
		SAMRecordIterator itr = blasrReader.iterator();
		while (itr.hasNext()) {
			blasrWriter.addAlignment(itr.next());
		}
		itr.close();
		blasrReader.close();
		blasrWriter.close();
		
	}
	
}
