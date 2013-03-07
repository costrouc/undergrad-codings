import java.io.Serializable;

/* 
 * Stores configuration data
 */
public class Configuration implements Serializable{
	private static final long serialVersionUID = 1L;
	private String emailAddress;
	private String ipAddress;
	
	public Configuration(String emailAddress,String ipAddress){
		this.emailAddress = emailAddress;
		this.ipAddress = ipAddress;
	}
	
	public Configuration(){
		this(null,null);
	}
	
	public void setEmailAddress(String emailAddress){
		this.emailAddress = emailAddress;
	}
	
	public void setipAddress(String ipAddress){
		this.ipAddress = ipAddress;
	}
	
	public String getemailAddress(){
		return emailAddress;
	}
	
	public String getipAddress(){
		return ipAddress;
	}
}
