import streamlit as st
import pandas as pd

# To upload a csv file
uploaded_file = st.sidebar.file_uploader(label="**Load your csv file**", 
                                         type = "csv", 
                                         label_visibility="hidden")

st.title("Nuclear cell step")
st.subheader("For cells with 2 components")
st.markdown("""---""")

# If a file has been uploaded
if uploaded_file is not None:

    # We read the file
    df = pd.read_csv(uploaded_file, sep=";")

    # We search the column "cell cycle" but if not used in the csv file, you can choose the correct one.
    if "cell cycle" in df.columns:
        cell_info = "cell cycle"
    else:
        cell_info = st.sidebar.selectbox(
            "Select the column that hold the information of the cell cycle of the nuclei",
            (df.columns),
            index=None,
            placeholder="Select cell cycle column (G1/G2/S)...",
            label_visibility='hidden')
    
    # Until a correct column is chosen
    if cell_info != None:

        # We group the cells on the Cell ID and ccocunt the number of nuclei
        groupby_count = df.groupby(["Parent Object ID (MO)"])[cell_info].count().reset_index()

        # We count the number of nuclei states by cells
        nuclei_present = list(groupby_count.value_counts(cell_info).index)
        # We count the number of cells with x nuclei
        cell_count_by_nuclei = list(groupby_count.value_counts(cell_info))
        # In a list we associate number of nuclei with number of cells
        counting_cell_nuclei = [[index, count] for index, count in zip(nuclei_present, cell_count_by_nuclei)]

        st.sidebar.markdown("""---""")

        # We give information about the distribution nucleis/cells
        for i in range(len(nuclei_present)):
            st.sidebar.write(f"""Number of nuclei 
                             {nuclei_present[i]}: 
                             {cell_count_by_nuclei[i]} cells
                             ({cell_count_by_nuclei[i]/sum(cell_count_by_nuclei):.1%})""")
        
        st.sidebar.markdown("""---""")
        
        # List containing all df with the different cells/nucleis
        list_df = []

        for n in nuclei_present:
            # We select the cells with n nuclei
            nuclei_id = list(groupby_count.loc[groupby_count[cell_info] == n, "Parent Object ID (MO)"])
            list_df.append(df[df["Parent Object ID (MO)"].isin(nuclei_id)])

        # User select which cells to look at
        nuclei_choice = st.selectbox(
            "Cells with N number of nuclei",
            (nuclei_present),
            index=None,
            placeholder="Select N nucleis in cells...")
        
        # To retrieve the dataframe from the list
        index = [i for i, n in enumerate(nuclei_present) if n == nuclei_choice][0]

        # We select the dataframe
        mask_df = list_df[index]
        mask_G1 = mask_df[mask_df[cell_info] == "G1"]
        mask_G2 = mask_df[mask_df[cell_info] == "G2"]
        mask_S = mask_df[mask_df[cell_info] == "S"]

        st.write("**Nuclei cell cycle**")
        st.write("Nuclei in G1: %.1f%% (%.1f%% total)" %(mask_G1.shape[0]/mask_df.shape[0]*100, mask_G1.shape[0]/df.shape[0]*100)) 
        st.write("Nuclei in G2: %.1f%% (%.1f%% total)" %(mask_G2.shape[0]/mask_df.shape[0]*100, mask_G2.shape[0]/df.shape[0]*100)) 
        st.write("Nuclei in S: %.1f%% (%.1f%% total)" %(mask_S.shape[0]/mask_df.shape[0]*100, mask_S.shape[0]/df.shape[0]*100)) 
        st.markdown("""---""")
        
        if nuclei_choice == 2:

            # If the number of nuclei is different than 2
            not_2_nuclear = list(groupby_count.loc[groupby_count[cell_info] != 2, "Parent Object ID (MO)"])

            df_2_nuclear = df[~df["Parent Object ID (MO)"].isin(not_2_nuclear)]

            df_not_2_nuclear = df[df["Parent Object ID (MO)"].isin(not_2_nuclear)]

            total_nuclei = mask_df.shape[0]

            similar_cell_cycle = mask_df[mask_df.duplicated(["Parent Object ID (MO)", cell_info], keep = False)]
            not_similar_cell_cycle = mask_df[~mask_df.duplicated(["Parent Object ID (MO)", cell_info], keep = False)]

            st.sidebar.write("Cells with similar nuclear cell cycle: %.1f%%" %(similar_cell_cycle.shape[0]/total_nuclei*100)) 
            st.sidebar.write("Cells with not similar nuclear cell cycle: %.1f%%" %(not_similar_cell_cycle.shape[0]/total_nuclei*100)) 
            st.sidebar.markdown("""---""")

            df_G1_G1 = similar_cell_cycle[similar_cell_cycle[cell_info] == "G1"]
            df_G2_G2 = similar_cell_cycle[similar_cell_cycle[cell_info] == "G2"]
            df_S_S = similar_cell_cycle[similar_cell_cycle[cell_info] == "S"]

            st.write("**Cells with similar nuclear cell cycle**")
            st.write("Cells with double G1/G1 nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G1_G1.shape[0]/total_nuclei*100, df_G1_G1.shape[0]/similar_cell_cycle.shape[0]*100)) 
            st.write("Cells with double G2/G2 nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G2_G2.shape[0]/total_nuclei*100, df_G2_G2.shape[0]/similar_cell_cycle.shape[0]*100)) 
            st.write("Cells with double S/S nuclear cell cycle: %.1f%% (%.1f%%)" %(df_S_S.shape[0]/total_nuclei*100, df_S_S.shape[0]/similar_cell_cycle.shape[0]*100)) 

            lists_dict = {"G1_G2": [], "G1_S": [], "G2_S": []}

            for id in not_similar_cell_cycle["Parent Object ID (MO)"].unique():
                mask = not_similar_cell_cycle[not_similar_cell_cycle["Parent Object ID (MO)"] == id] 
                list_c_cycle = list(mask[cell_info])
                
                if ("G1" in list_c_cycle and "G2" in list_c_cycle) or ("G2" in list_c_cycle and "G1" in list_c_cycle):
                    lists_dict["G1_G2"].append(id)
                elif ("G1" in list_c_cycle and "S" in list_c_cycle) or ("S" in list_c_cycle and "G1" in list_c_cycle):
                    lists_dict["G1_S"].append(id)
                elif ("G2" in list_c_cycle and "S" in list_c_cycle) or ("S" in list_c_cycle and "G2" in list_c_cycle):
                    lists_dict["G2_S"].append(id)

            df_G1_G2 = not_similar_cell_cycle[not_similar_cell_cycle["Parent Object ID (MO)"].isin(lists_dict["G1_G2"])]
            df_G1_S = not_similar_cell_cycle[not_similar_cell_cycle["Parent Object ID (MO)"].isin(lists_dict["G1_S"])]
            df_G2_S = not_similar_cell_cycle[not_similar_cell_cycle["Parent Object ID (MO)"].isin(lists_dict["G2_S"])]
            
            st.markdown("""---""")
            st.write("**Cells with not similar nuclear cell cycle**")
            st.write("Cells with G1/G2 nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G1_G2.shape[0]/total_nuclei*100, df_G1_G2.shape[0]/not_similar_cell_cycle.shape[0]*100)) 
            st.write("Cells with G1/S nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G1_S.shape[0]/total_nuclei*100, df_G1_S.shape[0]/not_similar_cell_cycle.shape[0]*100)) 
            st.write("Cells with G2/S nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G2_S.shape[0]/total_nuclei*100, df_G2_S.shape[0]/not_similar_cell_cycle.shape[0]*100)) 

            option = st.selectbox(
                "Which dataframe do you want to check?",
                ("Initial", "2 nuclears", "Not conform", "G1/G1", "G2/G2", "S/S", "G1/G2", "G1/S", "G2/S"),
                index=None,
                placeholder="Select a dataframe...")
            
            if option == "Initial":
                st.write(df)
            elif option == "2 nuclears":
                st.write(df_2_nuclear)
            elif option == "Not conform":
                st.write(df_not_2_nuclear)
            elif option == "G1/G1":
                st.write(df_G1_G1)
            elif option == "G2/G2":
                st.write(df_G2_G2)
            elif option == "S/S":
                st.write(df_S_S)
            elif option == "G1/G2":
                st.write(df_G1_G2)
            elif option == "G1/S":
                st.write(df_G1_S)
            elif option == "G2/S":
                st.write(df_G2_S)

    