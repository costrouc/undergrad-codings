import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;

/*
 * The main frame in the program allows for add,,edit, and remove to contact list
 * Allows edit all file options along with about
 */
public class MainFrame extends JFrame{
	private JMenuBar Window;
	private JFrame frame = this;
	private JPanel TablePanel = new JPanel();
	private JPanel ButtonPanel = new JPanel();
	private JTable Table;
	private JScrollPane scrollPane;
	private JMenu FileMenu, ConfigMenu, Help;
	private JMenuItem Exit, Configure, About;
	private JButton AddButton, EditButton, DeleteButton;

	private static final long serialVersionUID = 1L;

	public MainFrame(String Title) throws IOException{
		super(Title);
		SetUpFrame();
		SetUpMenuBar();	
		SetUpButtons();
		SetUpTable();
		
		getContentPane().add(TablePanel,"North");
		getContentPane().add(ButtonPanel,"South");
		
		setJMenuBar(Window);
		setVisible(true);
	}
	
	protected void SetUpFrame(){
		setLayout(new BorderLayout());
		setSize(800, 500);
		setResizable(false);
		ImageIcon image = new ImageIcon("emailicon.png");
		setIconImage(image.getImage());
	}
	
	protected void SetUpMenuBar(){
		Window = new JMenuBar();
		FileMenu = new JMenu("File");
	    Exit = new JMenuItem("Exit");
	    ConfigMenu = new JMenu("Configuration");
	    Configure = new JMenuItem("Configure");
	    Help = new JMenu("Help");
	    About = new JMenuItem("About");
	    Help.add(About);
	    addWindowListener(new WindowAdapter(){
	    	public void windowClosing(WindowEvent we){
				try {
					DataStore.getInstance().writeConfiguration();
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				try {
					DataStore.getInstance().writeContacts();
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	    	}
	    });
	    About.addActionListener(new ActionListener(){
	    	@Override
			public void actionPerformed(ActionEvent e) {
				new SystemInformationDlg(frame);
			}
	    });
	    	    
	    ConfigMenu.add(Configure);
	    Configure.addActionListener(new ActionListener(){
	    	@Override
	    	public void actionPerformed(ActionEvent e){
	    		Configuration config = DataStore.getInstance().getConfiguration();
	    		
	    		new ConfigurationDlg(frame, config);
	    	}
	    });    
	    
	    FileMenu.add(Exit);
	    Exit.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e) {
				if(e.getSource() == Exit){
					try {
						DataStore.getInstance().writeConfiguration();
					} catch (FileNotFoundException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					} catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					try {
						DataStore.getInstance().writeContacts();
					} catch (FileNotFoundException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					} catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					System.exit(0);
				}
			}
	    });
		
	    Window.add(FileMenu);
		Window.add(ConfigMenu);
		Window.add(Help);
	}
	
	protected void SetUpButtons(){
		AddButton = CreateButtonOnPanel("Add", 210, 400, 90, 30);
		AddButton.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e){
				Contact contact = new Contact();
				new ContactEditingDlg(frame,Table,contact);
				if (contact.getEmail().equals("") && 
						contact.getFirst().equals("") &&
						contact.getLast().equals("") &&
						contact.getPhone().equals("") &&
						contact.getPost().equals("")){
				}
				else{
					DataStore.getInstance().getContacts().addContact(contact);
				}
				EditButton.updateUI();
				DeleteButton.updateUI();
				AddButton.updateUI();
			}
		});
		
		EditButton = CreateButtonOnPanel("Edit", 340, 400, 90, 30);
		EditButton.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e){
				int i = Table.getSelectedRow();
				if (i != -1){
					Contact contact = DataStore.getInstance().getContacts().getContact(i);
					new ContactEditingDlg(frame,Table,contact);
					if (contact.getEmail().equals("") && 
							contact.getFirst().equals("") &&
							contact.getLast().equals("") &&
							contact.getPhone().equals("") &&
							contact.getPost().equals("")){
						DataStore.getInstance().getContacts().removeContact(i);
					}
					EditButton.updateUI();
					DeleteButton.updateUI();
					AddButton.updateUI();
				}
			}
		});
		
		DeleteButton = CreateButtonOnPanel("Delete", 470, 400, 90, 30);
		DeleteButton.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e){
				int i = Table.getSelectedRow();
				if (i != -1){
					 int reply = JOptionPane.showConfirmDialog(null, "Would you really like to delete?", "Delete Contact", JOptionPane.YES_NO_OPTION);
				        if (reply == JOptionPane.YES_OPTION) {
				        	DataStore.getInstance().getContacts().removeContact(i);
				        }
					
					Table.updateUI();
					EditButton.updateUI();
					DeleteButton.updateUI();
					AddButton.updateUI();
					Table.clearSelection();
				}
			}
		});
	}
	
	protected void SetUpTable(){
		Table = new JTable(new ContactTableModel());
		Table.getColumnModel().getColumn(0).setHeaderValue("First Name");
		Table.getColumnModel().getColumn(1).setHeaderValue("Last Name");
		Table.getColumnModel().getColumn(2).setHeaderValue("Postal Address");
		Table.getColumnModel().getColumn(3).setHeaderValue("Phone Number");
		Table.getColumnModel().getColumn(4).setHeaderValue("Email Address");
		scrollPane = new JScrollPane(Table, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		Table.setFillsViewportHeight(true);
		scrollPane.setPreferredSize(new Dimension(800,300));
		TablePanel.add(scrollPane);
		
        
		
		Table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		
		Table.addMouseListener(new MouseAdapter() {
			   public void mouseClicked(MouseEvent e) {
			      if (e.getClickCount() == 2) {
			         JTable target = (JTable) e.getSource();
			         int row = target.getSelectedRow();
			         if (row != -1){
			        	 Configuration config = DataStore.getInstance().getConfiguration();
			        	 Contact contact = DataStore.getInstance().getContacts().getContact(row);
			        	 new EmailTransmissionDlg(frame, config, contact.getEmail());
			         }
			         }
			   }
		});
	}
	
	protected JButton CreateButtonOnPanel(String name,int x1, int y1, int dx, int dy){
		JButton tempButton = new JButton(name);
		ButtonPanel.add(tempButton);
		tempButton.setBounds(x1,y1,dx,dy);
		return tempButton;
	}
}
