import streamlit as st
import pandas as pd

def get_column(df, expected_col, sidebar, placeholder):
    """
    Simplifies the selection of columns based on their presence or allowing the user to select.
    """
    if expected_col in df.columns:
        return expected_col
    else:
        return sidebar.selectbox(
            f"Select the column that holds the information of {placeholder}",
            df.columns,
            index=None,
            placeholder=f"Select {placeholder}...",
            label_visibility='hidden'
        )

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
    cell_info = get_column(df, "cell cycle", st.sidebar, "cell cycle of the nuclei")
    # User select the experiment column
    experiment_column = get_column(df, "samples", st.sidebar, "experiment")
    
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
            index=nuclei_present[0], # We put the first value
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
        # We retrieve the cell steps found in the dataset
        cell_steps = sorted(mask_df[cell_info].unique())

        #########################################################

        st.write("**Nuclei cell cycle distribution**")

        for cell_step in cell_steps:
            mask_cell_step = mask_df[mask_df[cell_info] == cell_step]

            try:
                # We write the % of cells in this condition in the cell step compare to total experiment
                st.write("Nuclei in %s: %.1f%% (%.1f%% total)" %(cell_step, mask_cell_step.shape[0]/mask_df.shape[0]*100, mask_cell_step.shape[0]/df.shape[0]*100)) 
            except ZeroDivisionError:
                # Handle the case where division by zero occurs
                st.error("No cells in this condition")

        #########################################################

        st.markdown("""---""")

        total_synchronized = 0

        # Create a DataFrame with unique "Parent Object ID (MO)"
        df_count = pd.DataFrame(mask_df["Parent Object ID (MO)"].unique(), columns=["Parent Object ID (MO)"])

        for cell_step in cell_steps:
            ### We count the occurence for each cell
            mask_cell_step = mask_df.loc[mask_df[cell_info] == cell_step, ["Parent Object ID (MO)", cell_info]]
            # Group by "Parent Object ID (MO)" and count occurrences of "G1", "G2", "S"
            cell_step_count = mask_cell_step.groupby("Parent Object ID (MO)").count().reset_index()
            
            df_count = df_count.merge(cell_step_count, how="left")
            
            rename_column = cell_step + "_count"

            # Rename the column
            df_count.rename(columns={cell_info: rename_column}, inplace=True)
            
            # Fill NaN values with 0
            df_count.fillna(0, inplace=True)
            
            # Convert count to integer type
            df_count[rename_column] = df_count[rename_column].astype(int)
            
            # To count the number of synchronized cells
            total_synchronized += df_count[df_count[rename_column] == nuclei_choice].shape[0]
        
        #########################################################

        st.write(f"**Cells with synchronized nucleis ({(total_synchronized/df_count.shape[0]):.1%})**")

        for cell_step in cell_steps:

            rename_column = cell_step + "_count"

            # We count the number of cells with totally synchronized nucleis
            sync_cell_step_count = df_count[df_count[rename_column] == nuclei_choice].shape[0]
            
            try:
                st.write("Nuclei in %s: %.1f%%" %(cell_step, sync_cell_step_count/df_count.shape[0]*100)) 
            except ZeroDivisionError:
                # Handle the case where division by zero occurs
                st.error("No cells in this condition")

        #########################################################

        st.markdown("""---""")

        st.write(f"**Cells with asynchronized nucleis ({1 - (total_synchronized/df_count.shape[0]):.1%})**")

        # Initialize the condition string
        condition = ""
        
        # Construct the condition string dynamically
        for cell_step in cell_steps:
            condition += f"(df_count['{cell_step + '_count'}'] != {nuclei_choice}) & "
        
        # Remove the trailing ' & ' from the condition string
        condition = condition[:-3]
        
        # Apply the condition and select the cell with asynchronous nucleis
        async_count_df = df_count[eval(condition)]

        for cell_step in cell_steps:

            # We calculate the proportion for each
            async_count_df.loc[:, [cell_step + "_count"]] = async_count_df.loc[:, [cell_step + "_count"]] / nuclei_choice
            # We measure the mean and standard deviation
            async_cell_step_mean = async_count_df[cell_step + "_count"].mean()
            async_cell_step_std = async_count_df[cell_step + "_count"].std()
            
            try:
                st.write(f"Proportion of nuclei in {cell_step}: {async_cell_step_mean:.1%} ± {async_cell_step_std:.1%}") 
            except ZeroDivisionError:
                # Handle the case where division by zero occurs
                st.error("No cells in this condition")

        #########################################################
                
        st.markdown("""---""")

        st.write(f"**Enrichment of cell cycle**")

        for cell_step in cell_steps:

            # We measure the mean and standard deviation
            enrich_cell_step_mean = (df_count/nuclei_choice)[cell_step + "_count"].mean()
            enrich_cell_step_std = (df_count/nuclei_choice)[cell_step + "_count"].std()
            
            try:
                st.write(f"Proportion of nuclei in G1: {enrich_cell_step_mean:.1%} ± {enrich_cell_step_std:.1%}") 
            except ZeroDivisionError:
                # Handle the case where division by zero occurs
                st.error("No cells in this condition")

        #########################################################
                
        st.markdown("""---""")


        
        

    