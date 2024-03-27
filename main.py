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

        # We group the cells on the Cell ID and count the number of nuclei
        groupby_count = df.groupby(["Parent Object ID (MO)", experiment_column])[cell_info].count().reset_index()

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

        # We give information about the distribution nucleis/experiments
        st.sidebar.subheader("Distribution of N nucleis by experiments")

        for exp in sorted(df[experiment_column].unique()):
            # We filter on the experiment and then measure mean and std
            mask = groupby_count[groupby_count.samples == exp]
            st.sidebar.write(f"""{exp}: 
                             {mask[cell_info].mean():.1f} 
                             ±{mask[cell_info].std():.1f}
                            ({(mask.shape[0] / df.shape[0]):.1%} cells)
                             """)

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
        
        # User select the experiment reference, e.g. control
        experiment_reference = st.selectbox(
            "Condition reference (e.g control):",
            (sorted(df[experiment_column].unique())),
            index=None,
            placeholder="Select experiments...")
        
        # User select the experiment condition he wants to see
        experiment_compare = st.selectbox(
            "Condition to compare with:",
            (sorted(df[experiment_column].unique())),
            index=None,
            placeholder="Select experiments...")
        
        st.markdown("""---""")
        
        # To retrieve the dataframe from the list
        index = [i for i, n in enumerate(nuclei_present) if n == nuclei_choice][0]

        # We select the reference dataframe
        mask_df_ref = list_df[index] if experiment_reference == None else list_df[index][list_df[index][experiment_column] == experiment_reference]
        mask_df_comp = None if experiment_compare == None else list_df[index][list_df[index][experiment_column] == experiment_compare]
        # We retrieve the cell steps found in the dataset
        cell_steps = sorted(mask_df_ref[cell_info].unique())

        #########################################################

        st.write("**Nuclei cell cycle distribution**")

        for cell_step in cell_steps:
            try:
                mask_cell_step_ref = mask_df_ref[mask_df_ref[cell_info] == cell_step]
                ref_percentage = mask_cell_step_ref.shape[0] / mask_df_ref.shape[0]
                total_ref_percentage = mask_cell_step_ref.shape[0]/df.shape[0]

                if isinstance(mask_df_comp, pd.DataFrame):
                    mask_cell_step_comp = mask_df_comp[mask_df_comp[cell_info] == cell_step]
                    # Calculate percentages
                    comp_percentage = mask_cell_step_comp.shape[0] / mask_df_comp.shape[0]
                    total_comp_percentage = mask_cell_step_comp.shape[0] / df.shape[0]
                    # Determine color based on comparison
                    color = "green" if comp_percentage >= ref_percentage else "red"
                    # We write the % of cells in this condition in the cell step compare to total experiment
                    
                    st.write(f"""Nuclei in {cell_step}: :{color}[{comp_percentage:.1%}] ({total_comp_percentage:.1%} total) 
                            VS {ref_percentage:.1%} ({total_ref_percentage:.1%} total)""")
                else:
                    st.write(f"""Nuclei in {cell_step}: {ref_percentage:.1%} ({total_ref_percentage:.1%} total)""")
            except ZeroDivisionError:
                # Handle the case where division by zero occurs
                st.error("No cells in this condition")

        #########################################################

        st.markdown("""---""")

        total_synchronized_ref = 0
        total_synchronized_comp = 0

        # Create a DataFrame with unique "Parent Object ID (MO)"
        df_count_ref = pd.DataFrame(mask_df_ref["Parent Object ID (MO)"].unique(), columns=["Parent Object ID (MO)"])

        for cell_step in cell_steps:
            ### We count the occurence for each cell
            mask_cell_step = mask_df_ref.loc[mask_df_ref[cell_info] == cell_step, ["Parent Object ID (MO)", cell_info]]
            # Group by "Parent Object ID (MO)" and count occurrences of "G1", "G2", "S"
            cell_step_count = mask_cell_step.groupby("Parent Object ID (MO)").count().reset_index()
            
            df_count_ref = df_count_ref.merge(cell_step_count, how="left")
            
            rename_column = cell_step + "_count"

            # Rename the column
            df_count_ref.rename(columns={cell_info: rename_column}, inplace=True)
            
            # Fill NaN values with 0
            df_count_ref.fillna(0, inplace=True)
            
            # Convert count to integer type
            df_count_ref[rename_column] = df_count_ref[rename_column].astype(int)
            
            # To count the number of synchronized cells
            total_synchronized_ref += df_count_ref[df_count_ref[rename_column] == nuclei_choice].shape[0]

        total_synchronized_ref_pourc = (total_synchronized_ref/df_count_ref.shape[0])

        if isinstance(mask_df_comp, pd.DataFrame):
            df_count_comp = pd.DataFrame(mask_df_comp["Parent Object ID (MO)"].unique(), columns=["Parent Object ID (MO)"])

            for cell_step in cell_steps:
                ### We count the occurence for each cell
                mask_cell_step = mask_df_comp.loc[mask_df_comp[cell_info] == cell_step, ["Parent Object ID (MO)", cell_info]]
                # Group by "Parent Object ID (MO)" and count occurrences of "G1", "G2", "S"
                cell_step_count = mask_cell_step.groupby("Parent Object ID (MO)").count().reset_index()
                
                df_count_comp = df_count_comp.merge(cell_step_count, how="left")
                
                rename_column = cell_step + "_count"

                # Rename the column
                df_count_comp.rename(columns={cell_info: rename_column}, inplace=True)
                
                # Fill NaN values with 0
                df_count_comp.fillna(0, inplace=True)
                
                # Convert count to integer type
                df_count_comp[rename_column] = df_count_comp[rename_column].astype(int)
                
                # To count the number of synchronized cells
                total_synchronized_comp += df_count_comp[df_count_comp[rename_column] == nuclei_choice].shape[0]

            total_synchronized_comp_pourc = (total_synchronized_comp/df_count_comp.shape[0])

            # Determine color based on comparison
            color = "green" if total_synchronized_comp_pourc >= total_synchronized_ref_pourc else "red"

            st.write(f"**Cells with synchronized nucleis (:{color}[{total_synchronized_comp_pourc:.1%}] VS {total_synchronized_ref_pourc:.1%})**")
        
        #########################################################
        else: 
            st.write(f"**Cells with synchronized nucleis ({total_synchronized_ref:.1%})**")
        
        for cell_step in cell_steps:

            try:
                rename_column = cell_step + "_count"

                # We count the number of cells with totally synchronized nucleis
                sync_cell_step_count_ref = df_count_ref[df_count_ref[rename_column] == nuclei_choice].shape[0]
                sync_cell_step_count_ref_pourc = sync_cell_step_count_ref/df_count_ref.shape[0]

                if isinstance(mask_df_comp, pd.DataFrame):
                    sync_cell_step_count_comp = df_count_comp[df_count_comp[rename_column] == nuclei_choice].shape[0]
                    sync_cell_step_count_comp_pourc = sync_cell_step_count_comp/df_count_comp.shape[0]

                    # Determine color based on comparison
                    color = "green" if sync_cell_step_count_comp_pourc >= sync_cell_step_count_ref_pourc else "red"

                    st.write(f"Nuclei in {cell_step}: :{color}[{sync_cell_step_count_comp_pourc:.1%}] VS {sync_cell_step_count_ref_pourc:.1%}")
                
                else: 
                    st.write(f"Nuclei in {cell_step}: {sync_cell_step_count_ref_pourc:.1%}")

            except ZeroDivisionError:
                # Handle the case where division by zero occurs
                st.error("No cells in this condition")

        #########################################################

        st.markdown("""---""")

        st.write(f"**Cells with asynchronized nucleis ({1 - (total_synchronized_ref/df_count_ref.shape[0]):.1%})**")

        # Initialize the condition string
        condition = ""
        
        # Construct the condition string dynamically
        for cell_step in cell_steps:
            condition += f"(df_count_ref['{cell_step + '_count'}'] != {nuclei_choice}) & "
        
        # Remove the trailing ' & ' from the condition string
        condition = condition[:-3]
        
        # Apply the condition and select the cell with asynchronous nucleis
        async_count_df = df_count_ref[eval(condition)]

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

        st.write(f"**Enrichment of cell cycle /// similar to first stat**")

        for cell_step in cell_steps:
            
            try:
                # We measure the mean and standard deviation
                enrich_cell_step_ref_mean = (df_count_ref/nuclei_choice)[cell_step + "_count"].mean()
                enrich_cell_step_ref_std = (df_count_ref/nuclei_choice)[cell_step + "_count"].std()

                if isinstance(mask_df_comp, pd.DataFrame):
                    enrich_cell_step_comp_mean = (df_count_comp/nuclei_choice)[cell_step + "_count"].mean()
                    enrich_cell_step_comp_std = (df_count_comp/nuclei_choice)[cell_step + "_count"].std()

                    # Determine color based on comparison
                    color = "green" if enrich_cell_step_comp_mean >= enrich_cell_step_ref_mean else "red"

                    st.write(f"""Proportion of nuclei in {cell_step}: :{color}[{enrich_cell_step_comp_mean:.1%}] ± :{color}[{enrich_cell_step_comp_std:.1%}]
                              VS {enrich_cell_step_ref_mean:.1%} ± {enrich_cell_step_ref_std:.1%}""") 
                else:
                    st.write(f"Proportion of nuclei in {cell_step}: {enrich_cell_step_ref_mean:.1%} ± {enrich_cell_step_ref_std:.1%}") 
            except ZeroDivisionError:
                # Handle the case where division by zero occurs
                st.error("No cells in this condition")

        #########################################################
                
        st.markdown("""---""")


        
        

    