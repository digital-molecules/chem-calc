
-----------------------------------------------------------------------------------
This app was written using PyCharm.

Make sure you install the 'rdkit' and 'streamlit' libraries using the 'pip install' command.
-----------------------------------------------------------------------------------
For the app to run independently from the IDE, you may need to install rdkit with conda and add it to requirements.txt along with streamlit
-----------------------------------------------------------------------------------
As long as you run main.py through PyCharm's GUI and it returns 0, you can run the app using the typical commands ('streamlit run main.py') on the terminal. If PyCharm shows errors but the app runs fine, this is due to lack of configuration on PyCharm's part that would allow it to properly test Streamlit code. You can ignore these errors as long as 0 is returned. If test_chem.py runs normally, there should be no other issues.
-----------------------------------------------------------------------------------
The app is usable as is, but the future goal is to add amine group calculations, improve the UI (because my brother called it ugly) and add visual descriptors for each molecule/SMILES string

#Have fun with SMILES strings!
