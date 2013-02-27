#ifndef ELEMENT_H
#define ELEMENT_H

enum element_info_types
{
  ELEMENT_NUMBER,
  ELEMENT_WEIGTH,
  ELEMENT_NAME,
  ELEMENT_SYMBOL
};

class Element
{
 public:
  Element(int atomicNumber);
  ~Element();

  double Get(info_types i);
  string Get(info_types i);
  
 private:
  ReadElemenetInfo(string elementInfo);
  
  double _name;
  string _symbol;
  int _number;
  double _weight;
};

class ElementFactory
{
 public:
  ElementFactory();
  ~ElementFactory();
  
  Element* GetElement(int i);
  Element* GetElement(string name);
  int GetNumberElements();
  
 private:
  Element _elements[];
  
  void ReadElemenetTable(); //If this method fails then Factory Should not work
};
