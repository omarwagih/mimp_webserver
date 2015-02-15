package models;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
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
	public Boolean incl_cent;
	public String ps_data;
	public Double p_thresh;
	public Double log_thresh;
	public String model_data;
	public String html;
	public List<String> mut_error;
	public List<String> ps_error;
	public Boolean has_error;
	
	public Scan(){
		job_id = "";
		fasta_data = "";
		mut_data = "";
		ps_data = "";
		incl_ps = "no";
		incl_cent = false;
		p_thresh = 0.5;
		log_thresh = 1.0;
		model_data = "hconf";
		mut_error = new ArrayList();
		ps_error = new ArrayList();
		has_error = false;
	}
	
	public String toString(){
		return "job_id = " + job_id + "\n" +
				"fasta_data = " + fasta_data + "\n"+
				"mut_data = " + mut_data + "\n"+
				"incl_ps = " + incl_ps + "\n"+
				"ps_data = " + ps_data + "\n"+
				"p_thresh = " + p_thresh + "\n"+
				"log_thresh = " + log_thresh + "\n";
	}
	
	public String modelDataName(){
		if(model_data.equals("hconf-fam")) return("_fam");
		if(model_data.equals("lconf")) return("_newman");
		return("");
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		File f = new File("public/R/mimp_web.R");
        Logger.info(f.exists()+"");
        
        String jobdir = "public/jobs/"+job_id;
        (new File(jobdir)).mkdirs();
        
		try {
			
			Map<String, String> fa = Help.readFasta(fasta_data);
			for(String line : mut_data.trim().split("\n")){
				if(line.trim().isEmpty()) continue;
				//System.out.println("LINE:" + line);
				
				String[] sp = line.split("\\s+");
				String gene = sp[0];
				//System.out.println("\tGENE:" + gene);
				String mut = sp[1];
				//System.out.println("\tMUT:" + mut);
				String seq = fa.get(gene);
				
				int mut_pos = Integer.parseInt(mut.substring(1, mut.length()-1));
				//System.out.println("\tmut_pos:" + Integer.toString(mut_pos));
				
				if(seq == null) continue;
				if(mut_pos > seq.length()) continue;
				
				char exp_aa = mut.charAt(0);
				char obs_aa = fa.get(gene).charAt(mut_pos - 1);
				
				if(exp_aa != obs_aa){
					mut_error.add(gene + " " + mut + ": expected '" + exp_aa + "', found '" + obs_aa +"'");
				}
			}
			
			if(!mut_error.isEmpty()){
				Application.writeJsonJobToFile(this);
				return;
			}
			
			if(incl_ps.equals("yes")){
				for(String line : ps_data.trim().split("\n")){
					String[] sp = line.split("\\s+");
					String gene = sp[0];
					int ps_pos = Integer.parseInt(sp[1]);
					String seq = fa.get(gene);
					if(null != seq){
						String obs_aa = seq.charAt(ps_pos - 1) + "";
						if(!obs_aa.matches("^[S|T|Y]$")){
							ps_error.add(gene + " " + ps_pos + ": expected 'S', 'T', or 'Y', found '" + obs_aa +"'");
						}
					}
					
				}
				
				if(!ps_error.isEmpty()){
					Application.writeJsonJobToFile(this);
					return;
				}
			}
			
			
			
			File fasta_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp");
			FileUtils.writeStringToFile(fasta_file, fasta_data+"\n");
	        
	        File mut_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp"); 
	        FileUtils.writeStringToFile(mut_file, mut_data+"\n");
	        
	        String cmd = "Rscript public/R/mimp_web.R ";
	        cmd = cmd + "--fasta " + fasta_file.getAbsolutePath() + " ";
	        cmd = cmd + "--mut " + mut_file.getAbsolutePath() + " ";
	        if(incl_ps.equals("yes")){
	        	File phos_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp"); 
	            FileUtils.writeStringToFile(phos_file, ps_data+"\n");
	        	cmd = cmd + "--phos " + phos_file.getAbsolutePath() + " ";
	        }
	        cmd = cmd + "--pthresh " + p_thresh + " ";
	        cmd = cmd + "--logthresh " + log_thresh + " ";
	        cmd = cmd + "--jobid " + job_id + " ";
	        cmd = cmd + "--mdata " + model_data + " ";
	        if(incl_cent) cmd = cmd + "--cent yes ";
	        
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
	    	//Help.writeFile(mut_error.toString(), "err.txt");
	    	Application.writeJsonJobToFile(this);
		} catch (Exception e) {
			this.has_error = true;
			Application.writeJsonJobToFile(this);
			// TODO Auto-generated catch block
			e.printStackTrace();
			
		} 
	}
}
