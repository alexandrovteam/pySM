
def scrape_from_kegg(kegg_id,base_url='http://www.kegg.jp/dbget-bin/www_bget?',base_suffix=''):
    import bs4 as BeautifulSoup

    url='{}{}{}'.format(base_url,kegg_id,base_suffix)
    if url.startswith('http'):
        #url='http://www.genome.jp/dbget-bin/www_bget?{}'.format(kegg_id)
        page = urllib2.urlopen(url).read()
    else:
        page = open(url)
        
    soup = BeautifulSoup.BeautifulSoup(page, "html.parser")
    alltables = soup.find_all("table",{"border":"0","cellpadding":"2"})
    table_data={}
    table_data['ID']='{}'.format(kegg_id)
    table_data['Mol weight']='0'
    table_data['Formula']=''
    table_data['Name']=''
    table_data['Pathway']=''
    if alltables == []:
        print 'table fail'
    else:
        for table in alltables:
            for row in table.find_all('tr'):
                row_header =row.find('th')
                row_text=row.find('td')
                if row_header is not None:
                    table_data[row_header.text]=row_text.text.encode('utf-8').strip() 
    return table_data['ID'],table_data['Name'], table_data['Formula'], table_data['Mol weight'],table_data['Pathway']

def download_from_kegg(kegg_id,output_dir):
    url='http://www.kegg.jp/dbget-bin/www_bget?{}'.format(kegg_id)
    #url='http://www.genome.jp/dbget-bin/www_bget?{}'.format(kegg_id)
    page = urllib2.urlopen(url).read()
    with open('{}/{}.html'.format(output_dir,kegg_id),'w') as f_out:
        f_out.write(page)
