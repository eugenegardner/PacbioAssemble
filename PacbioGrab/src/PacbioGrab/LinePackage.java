package PacbioGrab;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class LinePackage {

	private int numReads;
	private File fqFile;
	private Map<String, Integer> readNames = new HashMap<String, Integer>(); 
	
	public LinePackage (int num, File fqFile, Map<String, Integer> readNames) {
		
		numReads = num;
		this.fqFile = fqFile;	
		this.readNames = readNames;
	}
	
	public int getNumReads() {
		return numReads;
	}
	public File getFqFile() {
		return fqFile;
	}

	public Map<String, Integer> getReadNames() {
		return readNames;
	}

}
