import java.io.Serializable;
import java.util.Vector;

/*
 * Stores a vector of contacts
 */

@SuppressWarnings("serial")
public class ContactList implements Serializable{
	private Vector<Contact> contactList;
	
	public ContactList(Vector<Contact> contactList){
		this.contactList = contactList;
	}
	
	public ContactList(){
		this.contactList = new Vector<Contact>();
	}
	
	public void addContact(Contact contact){
		contactList.add(contact);
	}
	
	public void removeContact(int i){
		contactList.remove(i);
	}
	
	public Contact getContact(int i){
		return contactList.get(i);
	}
	
	public int size(){
		return contactList.size(); 
	}
}
