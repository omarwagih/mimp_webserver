package models;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.UUID;

import org.apache.commons.io.FileUtils;

import controllers.Application;
import controllers.Help;

import play.Logger;

public class Scan implements Runnable {

	public String job_id;
	public String fasta_data;
	public String mut_data;
	public String incl_ps;
	public String ps_data;
	public Integer a_min;
	public Integer a_max;
	public String html;
	
	public Scan(){
		job_id = "foo";
		fasta_data = "";
		mut_data = "";
		ps_data = "";
		a_min = 85;
		a_max = 15;
	}
	
	public String toString(){
		return "job_id = " + job_id + "\n" +
				"fasta_data = " + fasta_data + "\n"+
				"mut_data = " + mut_data + "\n"+
				"incl_ps = " + incl_ps + "\n"+
				"ps_data = " + ps_data + "\n"+
				"a_min = " + a_min + "\n"+
				"a_max = " + a_max + "\n";
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		File f = new File("public/R/mimp_web.R");
        Logger.info(f.exists()+"");
        
        String jobdir = "public/jobs/"+job_id;
        (new File(jobdir)).mkdirs();
        
		try {
			File fasta_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp");
			FileUtils.writeStringToFile(fasta_file, fasta_data);
	        
	        File mut_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp"); 
	        FileUtils.writeStringToFile(mut_file, mut_data);
	        
	        String cmd = "Rscript public/R/mimp_web.R ";
	        cmd = cmd + "--fasta " + fasta_file.getAbsolutePath() + " ";
	        cmd = cmd + "--mut " + mut_file.getAbsolutePath() + " ";
	        if(incl_ps.equals("yes")){
	        	File phos_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp"); 
	            FileUtils.writeStringToFile(phos_file, ps_data);
	        	cmd = cmd + "--mut " + phos_file.getAbsolutePath() + " ";
	        }
	        cmd = cmd + "--amin " + a_min + " ";
	        cmd = cmd + "--amax " + a_max + " ";
	        cmd = cmd + "--jobid " + job_id + " ";
	        
	        
	        
	        StringBuilder shellOutput = new StringBuilder();
	        Process p = Runtime.getRuntime().exec(cmd);

	        // Read in output returned from shell
	        // #########FOR DEBUG ONLY######
	        BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
	        String line;
	        while ((line = in.readLine()) != null) {
	            shellOutput.append(line + "\n");
	        }
	        
	        Logger.info(cmd);
	        //Runtime.getRuntime().exec(cmd);

	    	//html = Help.readFile(Application.jobPath(job_id) + "html.txt");
	    	
	    	Application.writeJsonJobToFile(this);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}
}
