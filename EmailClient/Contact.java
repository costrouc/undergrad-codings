

import java.io.Serializable;

/*
 * Stores the first, last, post, phone, email of the contact 
 * gives getter and setter methods
 */
public class Contact implements Serializable{
	private static final long serialVersionUID = 1L;
	private String FirstName;
	private String LastName;
	private String PostAddress;
	private String PhoneNum;
	private String EmailAddress;
	
	public Contact(String FirstName, String LastName, String PostAddress, String PhoneNum, String EmailAddress){
		this.FirstName = FirstName;
		this.LastName = LastName;
		this.PostAddress = PostAddress;
		this.PhoneNum = PhoneNum;
		this.EmailAddress = EmailAddress;
	}
	
	public Contact(){
		this.FirstName = new String();
		this.LastName = new String();
		this.PostAddress = new String();
		this.PhoneNum = new String();
		this.EmailAddress = new String();
	}
	
	public String getFirst(){
		return FirstName;
	}
	
	public String getLast(){
		return LastName;
	}
	
	public String getPost(){
		return PostAddress;
	}
	
	public String getPhone(){
		return PhoneNum;
	}
	
	public String getEmail(){
		return EmailAddress;
	}
	
	public void setFirst(String FirstName){
		this.FirstName = FirstName;
	}
	
	public void setLast(String LastName){
		this.LastName = LastName;
	}
	
	public void setPost(String PostAddress){
		this.PostAddress = PostAddress;
	}
	
	public void setPhone(String PhoneNum){
		this.PhoneNum = PhoneNum;
	}
	
	public void setEmail(String EmailAddress){
		this.EmailAddress = EmailAddress;
	}
}
