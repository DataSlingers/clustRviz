# -*- coding: utf-8 -*-
import scrapy
from random import randint


class InaugTextSpider(scrapy.Spider):
    name = "inaug_text"
    allowed_domains = ["http://www.presidency.ucsb.edu/"]
    start_urls = (
'http://www.presidency.ucsb.edu/ws/index.php?pid=25800',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25801',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25803',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25804',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25806',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25807',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25808',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25810',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25811',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25813',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25814',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25815',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25817',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25819',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25821',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25822',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25823',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25825',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25827',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25828',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25830',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25832',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25833',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25834',
'http://www.presidency.ucsb.edu/ws/index.php?pid=21804',
'http://www.presidency.ucsb.edu/ws/index.php?pid=15349',
'http://www.presidency.ucsb.edu/ws/index.php?pid=16022',
'http://www.presidency.ucsb.edu/ws/index.php?pid=16607',
'http://www.presidency.ucsb.edu/ws/index.php?pid=13282',
'http://www.presidency.ucsb.edu/ws/index.php?pid=10856',
'http://www.presidency.ucsb.edu/ws/index.php?pid=8032',
'http://www.presidency.ucsb.edu/ws/index.php?pid=1941',
'http://www.presidency.ucsb.edu/ws/index.php?pid=4141',
'http://www.presidency.ucsb.edu/ws/index.php?pid=6575',
'http://www.presidency.ucsb.edu/ws/index.php?pid=38688',
'http://www.presidency.ucsb.edu/ws/index.php?pid=16610',
'http://www.presidency.ucsb.edu/ws/index.php?pid=54183',
'http://www.presidency.ucsb.edu/ws/index.php?pid=25853',
'http://www.presidency.ucsb.edu/ws/index.php?pid=58745',
'http://www.presidency.ucsb.edu/ws/index.php?pid=102827',
'http://www.presidency.ucsb.edu/ws/index.php?pid=120000',
    )

    def parse(self, response):
        filename = response.xpath("//title/text()").extract()
        filename = filename[0] + str(randint(100,999))
        speech = response.xpath('//span[@class="displaytext"]').extract()
        print(filename)
        with open(filename,'wb') as f:
            f.write(speech[0].encode('utf8'))

