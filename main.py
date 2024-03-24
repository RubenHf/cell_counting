import streamlit as st
import pandas as pd

# To upload a csv file
uploaded_file = st.sidebar.file_uploader(label="**Load your csv file**", 
                                         type = "csv", 
                                         label_visibility="hidden")

st.title("Nuclear cell step")
st.subheader("Distribution of nuclei by cells")
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
        
    # User select the experiment column
    if "samples" in df.columns:
        experiment_column = "samples"
    else:
        experiment_column = st.sidebar.selectbox(
            "Select the column that hold the information of experiment",
            (df.columns),
            index=None,
            placeholder="Select samples column...")
    
    # Until a correct column is chosen
    if cell_info != None and experiment_column != None:

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
        st.sidebar.subheader("Distribution of cells by number of nucleis")
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
        
        # User select the experiment condition he wants to see
        experiment_choice = st.selectbox(
            "Condition:",
            (sorted(df[experiment_column].unique())),
            index=None,
            placeholder="Select experiments...")
        
        # To retrieve the dataframe from the list
        index = [i for i, n in enumerate(nuclei_present) if n == nuclei_choice][0]

        # We select the dataframe
        mask_df = list_df[index] if experiment_choice == None else list_df[index][list_df[index][experiment_column] == experiment_choice]
        mask_G1 = mask_df[mask_df[cell_info] == "G1"]
        mask_G2 = mask_df[mask_df[cell_info] == "G2"]
        mask_S = mask_df[mask_df[cell_info] == "S"]

        st.write("**Nuclei cell cycle distribution**")
        try:
            st.write("Nuclei in G1: %.1f%% (%.1f%% total)" %(mask_G1.shape[0]/mask_df.shape[0]*100, mask_G1.shape[0]/df.shape[0]*100)) 
            st.write("Nuclei in G2: %.1f%% (%.1f%% total)" %(mask_G2.shape[0]/mask_df.shape[0]*100, mask_G2.shape[0]/df.shape[0]*100)) 
            st.write("Nuclei in S: %.1f%% (%.1f%% total)" %(mask_S.shape[0]/mask_df.shape[0]*100, mask_S.shape[0]/df.shape[0]*100)) 
        except ZeroDivisionError:
            # Handle the case where division by zero occurs
            st.error("No cells in this condition")

        st.markdown("""---""")
        
        ### We count the occurence for each cell
        g1_mask_df = mask_df.loc[mask_df[cell_info] == "G1", ["Parent Object ID (MO)", cell_info]]
        g2_mask_df = mask_df.loc[mask_df[cell_info] == "G2", ["Parent Object ID (MO)", cell_info]]
        s_mask_df = mask_df.loc[mask_df[cell_info] == "S", ["Parent Object ID (MO)", cell_info]]

        # Group by "Parent Object ID (MO)" and count occurrences of "G1", "G2", "S"
        g1_counts = g1_mask_df.groupby("Parent Object ID (MO)").count().reset_index()
        g2_counts = g2_mask_df.groupby("Parent Object ID (MO)").count().reset_index()
        s_counts = s_mask_df.groupby("Parent Object ID (MO)").count().reset_index()

        # Create a DataFrame with unique "Parent Object ID (MO)"
        df_count = pd.DataFrame(mask_df["Parent Object ID (MO)"].unique(), columns=["Parent Object ID (MO)"])

        # Join the count of "G1", "G2", "S" occurrences
        for c_counts, c_counts_col in zip([g1_counts, g2_counts, s_counts], ["G1_count", "G2_count", "S_count"]):
            df_count = df_count.merge(c_counts, how="left")
            
            # Rename the column
            df_count.rename(columns={cell_info: c_counts_col}, inplace=True)
            
            # Fill NaN values with 0
            df_count.fillna(0, inplace=True)
            
            # Convert count to integer type
            df_count[c_counts_col] = df_count[c_counts_col].astype(int)

        # We count the number of cells with totally synchronized nucleis
        sync_G1_count = df_count[df_count["G1_count"] == nuclei_choice].shape[0]
        sync_G2_count = df_count[df_count["G2_count"] == nuclei_choice].shape[0]
        sync_S_count = df_count[df_count["S_count"] == nuclei_choice].shape[0]
        total_synchronized = sync_G1_count + sync_G2_count + sync_S_count

        st.write(f"**Cells with synchronized nucleis ({total_synchronized/df_count.shape[0]:.1%})**")
        
        try:
            st.write("Nuclei in G1: %.1f%%" %(sync_G1_count/df_count.shape[0]*100)) 
            st.write("Nuclei in G2: %.1f%%" %(sync_G2_count/df_count.shape[0]*100)) 
            st.write("Nuclei in S: %.1f%%" %(sync_S_count/df_count.shape[0]*100)) 
        except ZeroDivisionError:
            # Handle the case where division by zero occurs
            st.error("No cells in this condition")

        st.markdown("""---""")

        st.write(f"**Cells with asynchronized nucleis ({1 - (total_synchronized/df_count.shape[0]):.1%})**")

        # We select the cell with asynchronous nucleis
        async_count = df_count[(df_count["G1_count"] != nuclei_choice) & (df_count["G2_count"] != nuclei_choice) & (df_count["S_count"] != nuclei_choice)]
        # We calculate the proportion for each
        async_count.loc[:, ["G1_count", "G2_count", "S_count"]] = async_count.loc[:, ["G1_count", "G2_count", "S_count"]] / nuclei_choice
        # We measure the mean and standard deviation
        async_G1_mean = async_count["G1_count"].mean()
        async_G1_std = async_count["G1_count"].std()
        async_G2_mean = async_count["G2_count"].mean()
        async_G2_std = async_count["G2_count"].mean()
        async_S_mean = async_count["S_count"].std()
        async_S_std = async_count["S_count"].mean()
        
        try:
            st.write(f"Proportion of nuclei in G1: {async_G1_mean:.1%} ± {async_G1_std:.1%}") 
            st.write(f"Proportion of nuclei in G2: {async_G2_mean:.1%} ± {async_G2_std:.1%}") 
            st.write(f"Proportion of nuclei in S: {async_S_mean:.1%} ± {async_S_std:.1%}") 
        except ZeroDivisionError:
            # Handle the case where division by zero occurs
            st.error("No cells in this condition")

        st.markdown("""---""")

        st.write(f"**Enrichment of cell cycle**")

        # We measure the mean and standard deviation
        enrich_G1_mean = (df_count/nuclei_choice)["G1_count"].mean()
        enrich_G1_std = (df_count/nuclei_choice)["G1_count"].std()
        enrich_G2_mean = (df_count/nuclei_choice)["G2_count"].mean()
        enrich_G2_std = (df_count/nuclei_choice)["G2_count"].mean()
        enrich_S_mean = (df_count/nuclei_choice)["S_count"].std()
        enrich_S_std = (df_count/nuclei_choice)["S_count"].mean()
        
        try:
            st.write(f"Proportion of nuclei in G1: {enrich_G1_mean:.1%} ± {enrich_G1_std:.1%}") 
            st.write(f"Proportion of nuclei in G2: {enrich_G2_mean:.1%} ± {enrich_G2_std:.1%}") 
            st.write(f"Proportion of nuclei in S: {enrich_S_mean:.1%} ± {enrich_S_std:.1%}") 
        except ZeroDivisionError:
            # Handle the case where division by zero occurs
            st.error("No cells in this condition")

        st.markdown("""---""")


        
        

    