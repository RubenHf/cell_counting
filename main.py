import streamlit as st
import pandas as pd

# To upload a csv file
uploaded_file = st.sidebar.file_uploader(label="**Load your csv file**", 
                                         type = "csv", 
                                         label_visibility="hidden")

st.title("Nuclear cell step")
st.subheader("For cells with 2 components")
st.markdown("""---""")

if uploaded_file is not None:

    # We read the file
    df = pd.read_csv(uploaded_file, sep=";")

    groupby_count = df.groupby(["Parent Object ID (MO)"])["cell cycle"].count().reset_index()

    not_2_nuclear = list(groupby_count.loc[groupby_count["cell cycle"] != 2, "Parent Object ID (MO)"])

    df_2_nuclear = df[~df["Parent Object ID (MO)"].isin(not_2_nuclear)]

    df_not_2_nuclear = df[df["Parent Object ID (MO)"].isin(not_2_nuclear)]

    st.sidebar.markdown("""---""")

    st.sidebar.write("%d objects not conform out of %d total" 
                     %(df.shape[0] - df_2_nuclear.shape[0], df.shape[0]))
    
    st.sidebar.markdown("""---""")
    
    similar_cell_cycle = df_2_nuclear[df_2_nuclear.duplicated(["Parent Object ID (MO)", "cell cycle"], keep = False)]
    not_similar_cell_cycle = df_2_nuclear[~df_2_nuclear.duplicated(["Parent Object ID (MO)", "cell cycle"], keep = False)]

    st.sidebar.write("Cells with similar nuclear cell cycle: %.1f%%" %(similar_cell_cycle.shape[0]/df_2_nuclear.shape[0]*100)) 
    st.sidebar.write("Cells with not similar nuclear cell cycle: %.1f%%" %(not_similar_cell_cycle.shape[0]/df_2_nuclear.shape[0]*100)) 
    st.sidebar.markdown("""---""")

    df_G1_G1 = similar_cell_cycle[similar_cell_cycle["cell cycle"] == "G1"]
    df_G2_G2 = similar_cell_cycle[similar_cell_cycle["cell cycle"] == "G2"]
    df_S_S = similar_cell_cycle[similar_cell_cycle["cell cycle"] == "S"]

    st.markdown("""---""")
    st.write("**Cells with similar nuclear cell cycle**")
    st.write("Cells with double G1/G1 nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G1_G1.shape[0]/df_2_nuclear.shape[0]*100, df_G1_G1.shape[0]/similar_cell_cycle.shape[0]*100)) 
    st.write("Cells with double G2/G2 nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G2_G2.shape[0]/df_2_nuclear.shape[0]*100, df_G2_G2.shape[0]/similar_cell_cycle.shape[0]*100)) 
    st.write("Cells with double S/S nuclear cell cycle: %.1f%% (%.1f%%)" %(df_S_S.shape[0]/df_2_nuclear.shape[0]*100, df_S_S.shape[0]/similar_cell_cycle.shape[0]*100)) 

    lists_dict = {"G1_G2": [], "G1_S": [], "G2_S": []}

    for id in not_similar_cell_cycle["Parent Object ID (MO)"].unique():
        mask = not_similar_cell_cycle[not_similar_cell_cycle["Parent Object ID (MO)"] == id] 
        list_c_cycle = list(mask["cell cycle"])
        
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
    st.write("Cells with G1/G2 nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G1_G2.shape[0]/df_2_nuclear.shape[0]*100, df_G1_G2.shape[0]/not_similar_cell_cycle.shape[0]*100)) 
    st.write("Cells with G1/S nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G1_S.shape[0]/df_2_nuclear.shape[0]*100, df_G1_S.shape[0]/not_similar_cell_cycle.shape[0]*100)) 
    st.write("Cells with G2/S nuclear cell cycle: %.1f%% (%.1f%%)" %(df_G2_S.shape[0]/df_2_nuclear.shape[0]*100, df_G2_S.shape[0]/not_similar_cell_cycle.shape[0]*100)) 

    option = st.selectbox(
        "Which dataframe do you want to check?",
        ("Initial", "2 nuclears", "Not conform", "G1/G1", "G2/G2", "S/S", "G1/G2", "G1/S", "G2/S"),
        index=None,
        placeholder="Select a dataframe...",
        label_visibility="hidden")
    
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

    