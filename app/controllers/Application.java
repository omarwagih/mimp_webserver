package controllers;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;
import java.util.UUID;

import org.apache.commons.io.FileUtils;

import com.fasterxml.jackson.databind.JsonNode;

import models.Scan;
import play.*;
import play.mvc.*;
import play.data.Form;

import play.libs.Json;

import views.html.*;

public class Application extends Controller {
	final static Form<Scan> myForm = Form.form(Scan.class);
	
	public static Result wait_job(String jobid) throws InterruptedException {
		return ok(wait.render(jobid));
    }
	
	public static Result help() {
        return ok(help.render());
    }
	
	public static Result about() {
        return ok(about.render());
    }
	public static Result download() {
        return ok(download.render());
    }
	
    public static Result index() {
        return ok(index.render());
    }
    public static Result scan() {
    	Form<Scan> newform = myForm.fill(new Scan());
        return ok(scan.render(newform));
    }
    
    public static Result show_job(String jobid) throws IOException{
    	Scan job = getJobFromJsonFile(jobid);
    	if(job == null) return TODO;

    	job.html = Help.readFile(Application.jobPath(job.job_id) + "html.txt");
    	return ok(scandone.render(job));
    }

    public static Result submit() throws IOException{
    	// Get form
        Form<Scan> filledForm = myForm.bindFromRequest();
        Scan job = filledForm.get();

    	//job.job_id = 
        String jobid = "foo";
        jobid = UUID.randomUUID().toString();
    	job.job_id = jobid;
    	
        Logger.info(job.toString());
        
//        File f = new File("public/R/mimp_web.R");
//        Logger.info(f.exists()+"");
//        
//        File fasta_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp"); 
//        FileUtils.writeStringToFile(fasta_file, job.fasta_data);
//        
//        File mut_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp"); 
//        FileUtils.writeStringToFile(mut_file, job.mut_data);
//        
//        String cmd = "RScript public/R/mimp_web.R ";
//        cmd = cmd + "--fasta " + fasta_file.getAbsolutePath() + " ";
//        cmd = cmd + "--mut " + mut_file.getAbsolutePath() + " ";
//        if(job.incl_ps.equals("yes")){
//        	File phos_file = File.createTempFile(UUID.randomUUID().toString(), ".tmp"); 
//            FileUtils.writeStringToFile(phos_file, job.ps_data);
//        	cmd = cmd + "--mut " + phos_file.getAbsolutePath() + " ";
//        }
//        cmd = cmd + "--amin " + job.a_min + " ";
//        cmd = cmd + "--amax " + job.a_max + " ";;
//        
//        
//        StringBuilder shellOutput = new StringBuilder();
//        Process p = Runtime.getRuntime().exec(cmd);
//
//        // Read in output returned from shell
//        // #########FOR DEBUG ONLY######
//        BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
//        String line;
//        while ((line = in.readLine()) != null) {
//            shellOutput.append(line + "\n");
//        }
//        Logger.info(cmd);
//    	job.html = shellOutput.toString();
//    	writeJsonJobToFile(job);
        new Thread(job). start ( );
        
        return redirect("/wait/"+jobid);
    	
    }
    
    public static Scan getJobFromJsonFile(String jobid) {
        try {
            String jsonString = Help.readFile( jobPath(jobid) + jobid+".json" );
            Scan job = Json.fromJson(Json.parse(jsonString), Scan.class);

            return job;
        } catch (Exception e) {
            Logger.error(e.getMessage());
            return null;
        }
    }

    public static boolean writeJsonJobToFile(Scan job) {
        try {
            JsonNode jn = Json.toJson(job);
            String path = jobPath(job.job_id) +job.job_id+".json";
            Logger.info("Writing job to: "+path);
            Help.writeFile(jn.toString(), path);
            return true;
        } catch (Exception e) {
            return false;
        }
    }
    
    public static String jobPath(String jobid){
    	return "public/jobs/" +jobid + "/";
    }
}
