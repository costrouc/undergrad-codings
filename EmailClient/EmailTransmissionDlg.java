import java.awt.GridLayout;
import javax.mail.*;
import javax.mail.internet.*;
import javax.mail.internet.MimeMessage.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;

// Class for Transmitting an Email using java mail. 
// Provides an editable body, subject, and to field
// Source field is only editable through system configuration
// Destination mail address is also verified
public class EmailTransmissionDlg extends JDialog implements ActionListener{

	private static final long serialVersionUID = 1L;
    private Configuration config;
	private JPanel emailPanel;
	
	private JLabel Source_label;
	private JLabel Email_label = new JLabel("To: ");
	private JLabel Subject_label = new JLabel(" Subject:");
	private JLabel Body_label = new JLabel("Body:    ");
	private JTextField Email_tfield = new JTextField(40);
	private JTextField Subject_tfield = new JTextField(40);
	private JTextArea Body_tfield = new JTextArea(40,10);
	private JButton send = new JButton("Send");
	private JButton exit = new JButton("Cancel");
	
	public EmailTransmissionDlg(JFrame f, Configuration config, String contact_mail){
		super(f);
		setModal(true);
		this.config = config;
		emailPanel = new JPanel();
		emailPanel.setLayout(new GridLayout(9,1));
		Body_tfield.setLineWrap(true);
		Body_tfield.setWrapStyleWord(true);
		if(config.getemailAddress() == null)
		    Source_label = new JLabel("Source Address: ");
		else
			 Source_label = new JLabel("Source Address: " + config.getemailAddress());
		emailPanel.add(Source_label);
		emailPanel.add(Email_label);
		Email_tfield.setText(contact_mail);
		emailPanel.add(Email_tfield);
		emailPanel.add(Subject_label);
		emailPanel.add(Subject_tfield);
		emailPanel.add(Body_label);
		JScrollPane scrollArea = new JScrollPane(Body_tfield);
		emailPanel.add(scrollArea);
		exit.addActionListener(this);
		send.addActionListener(this);
		emailPanel.add(send);
		emailPanel.add(exit);
		
		setSize(500,700);
		
		setContentPane(emailPanel);
		
		setVisible(true);
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		if(arg0.getSource() == exit)
		   this.dispose();  //disposes object if exit is clicked
		else if(arg0.getSource() == send){
			String ToField = Email_tfield.getText();
			StringTokenizer tokens = new StringTokenizer(ToField,",; ");
			while (tokens.hasMoreElements()){
			
			
			
			   // checks for a valid email address
			   String email = tokens.nextToken();
			   Pattern p = Pattern.compile(".+@.+\\.[a-z]+");
			   Matcher m = p.matcher(email);
			   boolean matchFound = m.matches();

			   if(!matchFound){
			      JOptionPane.showMessageDialog(null, "Invalid Email Address: " + email);
			   }else{
			       try {
				      Properties props = new Properties();
				      props.put("mail.smtp.host", config.getipAddress());
				
			          Session session = Session.getDefaultInstance(props, null);
				
				      Message msg = new MimeMessage(session); 
				      msg.setRecipient(RecipientType.TO, new InternetAddress(email));
				      msg.setSubject(Subject_tfield.getText());
				      msg.setFrom(new InternetAddress(config.getemailAddress()));
				      msg.setText(Body_tfield.getText());
				      Transport.send(msg);
				      System.out.println("Message sent");
			      } catch (Exception exc) {
				      System.out.println("Exception: " + exc);
			      }
			   }
			      this.dispose();
		   }
	   }
	}
}


