# -*- coding: utf-8 -*-
import scrapy
from random import randint


class ConvTextSpider(scrapy.Spider):
    name = "conv_text"
    allowed_domains = ["http://www.presidency.ucsb.edu"]
    start_urls = (
'http://www.presidency.ucsb.edu/ws/index.php?pid=117935',
    )

    def parse(self, response):
        filename = response.xpath("//title/text()").extract()
        filename = filename[0] + str(randint(100,999))
        speech = response.xpath('//span[@class="displaytext"]').extract()
        print(filename)
        with open(filename,'wb') as f:
            f.write(speech[0].encode('utf8'))

