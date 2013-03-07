import java.io.IOException;
// Created by Paul Bowen and Christopher Ostrouchov
// Date: 12/7/2012
// Generics were not used as we did not have to cast at any point
public class MainDriver {

	public static void main(String[] args) throws ClassNotFoundException, IOException {
		DataStore.getInstance().loadConfiguration();
		DataStore.getInstance().loadContacts();
		new MainFrame("Simple Mail");

	}

}
