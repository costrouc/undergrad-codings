import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.JTextField;


/*
 * In this dialog we allow a user to either edit or add a new contact
 * the fields are checked for validity if necissary
 */
public class ContactEditingDlg extends JDialog implements ActionListener{
	private static final long serialVersionUID = 1L;

	private Contact contactInfo;
	
	private JPanel contactPanel = new JPanel();
	
	private JTable Table;
	
	private JLabel firstNameLabel = new JLabel("First Name:");
	private JLabel lastNameLabel = new JLabel("Last Name:");
	private JLabel postalAddressLabel = new JLabel("Postal Address:");
	private JLabel phoneNumberLabel = new JLabel("Phone Number:");
	private JLabel emailAddressLabel = new JLabel("Email Address:");
	
	private JTextField firstNameTField;
	private JTextField lastNameTField;
	private JTextField postalAddressTField;
	private JTextField phoneNumberTField;
	private JTextField emailAddressTField;
	
	private JButton cancel;
	private JButton save;
	
	public ContactEditingDlg(JFrame f, JTable table, Contact contactInfo){
		super(f);
		setModal(true);
		this.contactInfo = contactInfo;
		this.Table = table;
		
		SetUpTextFields();
		SetUpButtons();
		SetUpPanel();
		
		
		setSize(600,200);
		setResizable(false);
		setContentPane(contactPanel);
		setVisible(true);
	}
	
	public void SetUpPanel(){
		contactPanel.setLayout(new GridLayout(6,2));
		
		contactPanel.add(firstNameLabel);
		contactPanel.add(firstNameTField);
		contactPanel.add(lastNameLabel);
		contactPanel.add(lastNameTField);
		contactPanel.add(postalAddressLabel);
		contactPanel.add(postalAddressTField);
		contactPanel.add(phoneNumberLabel);
		contactPanel.add(phoneNumberTField);
		contactPanel.add(emailAddressLabel);
		contactPanel.add(emailAddressTField);

		contactPanel.add(save);
		contactPanel.add(cancel);
	}
	
	public void SetUpTextFields(){
		firstNameTField = new JTextField();
		firstNameTField.setText(contactInfo.getFirst());
		lastNameTField = new JTextField();
		lastNameTField.setText(contactInfo.getLast());
		postalAddressTField = new JTextField();
		postalAddressTField.setText(contactInfo.getPost());
		phoneNumberTField = new JTextField();
		phoneNumberTField.setText(contactInfo.getPhone());
		emailAddressTField = new JTextField();
		emailAddressTField.setText(contactInfo.getEmail());
	}
	
	public void SetUpButtons(){
		cancel = new JButton("cancel");
		save = new JButton("save");
		save.addActionListener(this);
		cancel.addActionListener(this);
	}
	
	public void actionPerformed(ActionEvent e) {
		int test = 0;
		if (e.getSource() == cancel){
			this.dispose();
		}
		if (e.getSource() == save){
			//check validity.....
			
				contactInfo.setFirst(firstNameTField.getText());
				contactInfo.setLast(lastNameTField.getText());
				contactInfo.setPost(postalAddressTField.getText());
				contactInfo.setPhone(phoneNumberTField.getText());
				
				
				//Checks if this is a valid email address
				String email = emailAddressTField.getText();
				Pattern p = Pattern.compile(".+@.+\\.[a-z]+");
				Matcher m = p.matcher(email);
				boolean matchFound = m.matches();

				if(!matchFound){
					JOptionPane.showMessageDialog(null, "Invalid Email Address Entered");
					test = 1;
				}
				else{
					contactInfo.setEmail(emailAddressTField.getText());
				}

				Table.updateUI();
				if(test == 0)
				   this.dispose();
			}
		
	}
}
