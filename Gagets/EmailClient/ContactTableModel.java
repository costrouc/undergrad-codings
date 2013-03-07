import javax.swing.table.AbstractTableModel;

/*
 * Allows a model for the table in such that the JTable can be updated
 */

public class ContactTableModel extends AbstractTableModel{
	private static final long serialVersionUID = 1L;
	
	private int columnCount;
	private ContactList contactList;
	
	public ContactTableModel(){
		DataStore data = DataStore.getInstance();
		columnCount = 5;
		contactList = data.getContacts();
	}
	
	@Override
	public int getColumnCount() {
		return columnCount;
	}

	@Override
	public int getRowCount() {
		return contactList.size();
	}
	
	@Override
	public Object getValueAt(int row, int col) {
		switch (col){
		case 0:
			return contactList.getContact(row).getFirst();
		case 1:
			return contactList.getContact(row).getLast();
		case 2:
			return contactList.getContact(row).getPost();
		case 3:
			return contactList.getContact(row).getPhone();
		case 4:
			return contactList.getContact(row).getEmail();
		default:
			return new String("Bad Call");
		}
	}
}
