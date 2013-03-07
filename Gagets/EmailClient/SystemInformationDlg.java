import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

/*
 * About information
 */
public class SystemInformationDlg extends JDialog implements ActionListener{

	private static final long serialVersionUID = 1L;

	private JPanel aboutPanel;
	
	private JLabel authors = new JLabel("This software was developed by: Christopher Ostrouchov, Paul Bowen");
	private JLabel dateCreated = new JLabel("(c) November 29 2012");
	private JLabel about = new JLabel("This software was created to enable a lightweight email client");
	
	private JButton ok = new JButton("OK");
	
	public SystemInformationDlg(JFrame f){
		super(f);
		setModal(true);
		aboutPanel = new JPanel();
		aboutPanel.add(authors);
		aboutPanel.add(dateCreated);
		aboutPanel.add(about);
		ok.addActionListener(this);
		aboutPanel.add(ok);
		
		setSize(600,200);
		
		setContentPane(aboutPanel);
		
		setVisible(true);
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		this.dispose();	
	}
}
