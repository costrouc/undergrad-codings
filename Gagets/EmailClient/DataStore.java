import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;


/*
 * Singleton that maintains reads and writes to the disk as well at storing the config data and contact data
 */
public class DataStore {
	private static DataStore m_singleton = new DataStore();
	private ContactList contactList; 
	private Configuration EmailConfig;
	
	private DataStore(){
		contactList =  new ContactList();
		EmailConfig = new Configuration();
	}
	
	public static DataStore getInstance(){
		return m_singleton;
	}
	
	public void writeContacts() throws FileNotFoundException, IOException{
		ObjectOutputStream output = new ObjectOutputStream(new FileOutputStream("Contacts.txt"));
		output.writeObject(contactList);
		output.close();
	}
	
	public void loadContacts() throws IOException, ClassNotFoundException{
		if(!new File("Contacts.txt").exists()){
			return;
		}
		ObjectInputStream input = new ObjectInputStream(new FileInputStream("Contacts.txt"));
		contactList = (ContactList) input.readObject();
		input.close();
	}
	
	public void writeConfiguration() throws FileNotFoundException, IOException{
		ObjectOutputStream output = new ObjectOutputStream(new FileOutputStream("Configuration.txt"));
		output.writeObject(EmailConfig);
		output.close();
	}
	
	public void loadConfiguration() throws IOException, ClassNotFoundException{
		if(!new File("Configuration.txt").exists()){
			return;
		}
		ObjectInputStream input = new ObjectInputStream(new FileInputStream("Configuration.txt"));
		EmailConfig = (Configuration)input.readObject();
		input.close();
	}
	
	public ContactList getContacts(){
		return contactList;
	}
	
	public Configuration getConfiguration(){
		return EmailConfig;
	}
		
}
