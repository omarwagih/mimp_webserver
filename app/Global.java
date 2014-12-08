

import play.*;
import play.mvc.*;
import play.mvc.Http.*;
import play.libs.F.*;
import static play.mvc.Results.*;
import views.html.*;

public class Global extends GlobalSettings {
	
	public Promise<SimpleResult> onHandlerNotFound(RequestHeader request) {
        return Promise.<SimpleResult>pure(notFound(
        		views.html.oops.render("handlerNotFound")
        ));
    }
	
	 public Promise<SimpleResult> onBadRequest(RequestHeader request, String error) {
		 Logger.error(error);
        return Promise.<SimpleResult>pure(badRequest(
        		views.html.oops.render("badRequest")
        ));
	 }
	 
	 
	 public Promise<SimpleResult> onError(RequestHeader request, String error) {
		 Logger.error(error);
	        return Promise.<SimpleResult>pure(badRequest(
	        		views.html.oops.render("error")
	        ));
		 }
	 
	 public Promise<SimpleResult> onExecutionException(RequestHeader request, String error) {
		 	Logger.error(error);
	        return Promise.<SimpleResult>pure(badRequest(
	        		views.html.oops.render("error")
	        ));
		 }
	
}
