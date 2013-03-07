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
import javax.swing.JTextField;

/*
 * This dialog allows you to configure the stmp address and from email address
 * The validity of the configuration is also checked.
 */
public class ConfigurationDlg extends JDialog implements ActionListener{

	private static final long serialVersionUID = 1L;
	
	private Configuration config;
	
	private JPanel configPanel;
	
	private JLabel Email_label = new JLabel("Email Address: ");
	private JLabel SMTP_label = new JLabel("SMTP Server:    ");
	private JTextField email_tfield;
	private JTextField smtp_tfield;
	private JButton save = new JButton("Save");
	private JButton exit = new JButton("Cancel");
	
	public ConfigurationDlg(JFrame f, Configuration c){
		super(f);
		setModal(true);
		config = c;
		smtp_tfield = new JTextField(config.getipAddress(),30);
		email_tfield = new JTextField(config.getemailAddress(),30);
		configPanel = new JPanel();
		configPanel.add(Email_label);
		configPanel.add(email_tfield);
		configPanel.add(SMTP_label);
		configPanel.add(smtp_tfield);
		exit.addActionListener(this);
		save.addActionListener(this);
		configPanel.add(save);
		configPanel.add(exit);
		
		setSize(500,150);
		
		setContentPane(configPanel);
		
		setVisible(true);
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		int test = 0;
		if(arg0.getSource() == exit)
		   this.dispose();
		else if(arg0.getSource() == save){
			if(smtp_tfield.getText().contains("smtp."))
			   config.setipAddress(smtp_tfield.getText());
			else{
				JOptionPane.showMessageDialog(null, "Invalid SMTP Server Entered");
				test = 1;
			}
			String email = email_tfield.getText();
			Pattern p = Pattern.compile(".+@.+\\.[a-z]+");
			Matcher m = p.matcher(email);
			boolean matchFound = m.matches();

			if(!matchFound){
			   JOptionPane.showMessageDialog(null, "Invalid Email Address Entered");
			   test = 1;
			}else{
			config.setEmailAddress(email_tfield.getText());
			if(test == 0)
			   this.dispose();
			}
		}
	}
}



