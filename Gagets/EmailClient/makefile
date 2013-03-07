JFLAGS = -g
JC = javac
LIB = ".:mail.jar"
.SUFFIXES: .java .class
.java.class:
	$(JC) $(JFLAGS) -classpath $(LIB) $*.java 

CLASSES = ConfigurationDlg.java \
	Configuration.java \
	ContactEditingDlg.java \
	Contact.java \
	ContactList.java \
	ContactTableModel.java \
	DataStore.java \
	EmailTransmissionDlg.java \
	MainDriver.java \
	MainFrame.java \
	SystemInformationDlg.java

default: classes

classes: $(CLASSES:.java=.class)

clean:
	$(RM) *.class

test:
	java MainDriver