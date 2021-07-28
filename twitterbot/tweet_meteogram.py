from TwitterAPI import TwitterAPI
import json

def tweet_meteogram(pngfilepath,text):

    tokens = json.load(open("tokens.json"))

    api = TwitterAPI(tokens.get("api_key"),tokens.get("api_secret"),tokens.get("app_token"),tokens.get("app_secret"))

    pngfile = open(pngfilepath, 'rb')
    pngdata = pngfile.read()
    r = api.request('statuses/update_with_media', {'status':text}, {'media[]':pngdata})
    return r.status_code