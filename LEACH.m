clear;

%Alan boyutlarý x ve y metre olarak
xm = 100;
ym = 100;

%Baz istasyonun x ve y kordinatlarý
sink.x = 1.5 * xm;
sink.y = 0.5 * ym;

%Alandaki node sayýsý
n = 200;

%Bir nodun clusterhead seçilmesi için gereken olasýlýk
p = 0.2;

%Baþlangýç enerjileri Joule cinsinden
Eo = 0.1; %Her bir nodun baþlangýç enerjisi
%Eelec=Etx=Erx Veri iletmek veya almak için elektronik cihaz tarafýndan harcanan enerji
ETX = 50 * 0.000000001; %Veri iletim enerjisi
ERX = 50 * 0.000000001; %Veri alma enerjisi

Efs = 10 * 0.000000000001; %Daðýtým sýrasýnda kaybolan enerji enerjisi
Emp = 0.0013 * 0.000000000001; %Enerji transferinin amplifikasyon katsayýsý
EDA = 5 * 0.000000001; %CH tarafýndan veri toplama sýrasýnda harcadýklarý enerjisi

k = 4000; %Veri paket boyutu

%Advanced node yüzdesi(yüksek enerjili node yüzdesi)
m = 0.0;
%\alpha - advanced nodelara enerji atarken kullanýlan kat sayý
a = 1;

%Maksimum iterasyon sayýsý
rmax = 100;

%BS'den alýcý düðüme olan uzaklýk, Emp veya Efs'e baðlýdýr.
do = sqrt(Efs / Emp);

%Sensor networkün oluþturulmasý
figure(1);
hold off;
for i = 1:1:n
    S(i).xd = rand(1, 1) * xm; %i ninci nodun x kordinatýný  xm degeri aralýgýnda üretip (S node dizisi)
    XR(i) = S(i).xd; %XR dizisine x kordinatýný atadýk
    S(i).yd = rand(1, 1) * ym; %i ninci nodun y kordinatýný  yd degeri aralýgýnda üretip (S node dizisi)
    YR(i) = S(i).yd; %YR dizisine y kordinatýný atadýk
    S(i).G = 0; %Clusterhead (CH) olabilecek advencad(enerjisi en yüksek olanlar) node larýn yada normal nodelarýn kümesi
    S(i).type = 'N'; %Baþlangýçta hiç cluster head olmadý için nodelarýn tipini n olarak iþaretledik yani normal node
 
    %Belirtilen advance node oranýna göre, Normal node larýn seçilmesi ve enerjilerinin atanmasý
    if (i >= m * n + 1)
        S(i).E = Eo;
        S(i).IS_ADVANCED = 0; %Advanced olmayacak
        plot(S(i).xd, S(i).yd, 'o'); %Normal node u ekrana çizdirdik
        hold on;
    end
    %Belirtilen advance node oranýna göre ilk x tanesini advence node olarak seç ve enerjilerinin atanmasý
    if (i < m * n + 1)
        S(i).E = Eo * (1 + a)
        S(i).IS_ADVANCED = 1; %Advanced olacak
        plot(S(i).xd, S(i).yd, '+'); %Yüksek enerjili(advanced) nodu ekrana çizdirdik
        hold on;
    end
end

%S dizinde node listemiz var node dizisinin en sonuna  baz istesyonunun kordinatlarýný ekledik
S(n + 1).xd = sink.x;
S(n + 1).yd = sink.y;
plot(S(n + 1).xd, S(n + 1).yd, 'x'); %Baz istasyonunu ekrana x olarak çizdirdik

%Ýlk iterasyon için tanýmlar

%Cluster head sayacý
countCHs = 0;
%Round daki cluster head sayýsýný tutan degiþken
rcountCHs = 0;

%Ýlk ölümün gerçekleþip gerçekleþmediðinden haberdar oldugumuz degiþken
flag_first_dead = 0;

for r = 0:1:rmax
    %r
 
    %Operation for epoch
    if (mod(r, round(1 / p)) == 0)
        for i = 1:1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end
 
    hold off;
 
    %Her round baþladýðýnda toplam ölü node sayýsý tutuyor, her round baþladýðýnda sýfýrlanýyor
    total_dead = 0;
    %Her round baþladýðýnda Advanced Ölü node sayýsý, her round baþladýðýnda sýfýrlanýyor
    dead_a = 0;
    %Her round baþladýðýnda Normal ölü node sayýsý, her round baþladýðýnda sýfýrlanýyor
    dead_n = 0;
 
    %BS ve CH ye gönderilen paket sayacý
    packets_TO_BS = 0;
    packets_TO_CH = 0;
    %Roundlarda CH ve BS a gönderilen paket sayacý dizi
    PACKETS_TO_CH(r + 1) = 0;
    PACKETS_TO_BS(r + 1) = 0;
 
    figure(1);
 
    for i = 1:1:n
        %Ölü node var mý kontrol et ve say
        if (S(i).E <= 0)
            plot(S(i).xd, S(i).yd, 'red .'); %En son canlý kalanlar grafiginde ölü nodelarý kýrmýzý nokta olarak göster
            total_dead = total_dead + 1;
            if (S(i).IS_ADVANCED == 1)
                dead_a = dead_a + 1;
            end
            if (S(i).IS_ADVANCED == 0)
                dead_n = dead_n + 1;
            end
            hold on;
        end
     
        figure(2);
        %En son canlý kalan node lar
        if S(i).E > 0
            S(i).type = 'N';
            if (S(i).IS_ADVANCED == 0)
                plot(S(i).xd, S(i).yd, 'o');
            end
            if (S(i).IS_ADVANCED == 1)
                plot(S(i).xd, S(i).yd, '+');
            end
            hold on;
        end
    end
    plot(S(n + 1).xd, S(n + 1).yd, 'x'); %En son canlý kalanlar grafiginde baz istasyonunuda göster
 
    %STATISTICS dizisinin içerisine toplam ölü advanced ölü normal ölü sayýsýný her bir round için atadýk
    STATISTICS(r + 1).DEAD = total_dead;
    STATISTICS(r + 1).DEAD_ADVANCED = dead_a;
    STATISTICS(r + 1).DEAD_NORMAL = dead_n;
 
    %DEAD dizisinin içerisine toplam ölü advanced ölü normal ölü sayýsýný her bir round için atadýk
    DEAD(r + 1) = total_dead;
    DEAD_N(r + 1) = dead_n;
    DEAD_A(r + 1) = dead_a;
 
    %Ýlk nodun öldüðü roundu tespit ettik
    if (total_dead == 1)
        if (flag_first_dead == 0)
            first_dead = r;
            flag_first_dead = 1;
        end
    end
    %Cluster head sayacý, her round baþladýgýnda sýfýrlanýr
    countCHs = 0;
    cluster = 1; %Cluster sayýsýný 1 olarak atýyo ilk baþta
    for i = 1:1:n
        if (S(i).E > 0) %Elimizdeki node ölü degilse
            temp_rand = rand;
            if ((S(i).G) <= 0)
             
                %Clusterhead seçim
                if (temp_rand <= (p / (1 - p * mod(r, round(1 / p))))) %Clusterheadin seçim formülü
                    countCHs = countCHs + 1; %Cluster sayýsýný 1 artýr
                    packets_TO_BS = packets_TO_BS + 1; %Baz istasyonuna gönderilen paket sayýsýný 1 artýr
                    PACKETS_TO_BS(r + 1) = packets_TO_BS; %O round içerisinde BS ye gönderilen paket sayýsýný güncelliyor
                 
                    S(i).type = 'C'; %Seçilen nodun tipini cluster olarak ayarla
                    S(i).G = round(1 / p) - 1; %Cluster head olarak seçilen nodeun G deðerini ata
                    C(cluster).xd = S(i).xd; %Clusterýn x kordinatýný C dizisine ata
                    C(cluster).yd = S(i).yd; %Clusterýn y kordinatýný C dizisine ata
                    plot(S(i).xd, S(i).yd, 'k*'); %En son yani figure 10 a clusterheadlari * olarak çiz
                 
                    %Clusterhead olarak seçilen node ile baz istasyonunun birbirine olan uzaklýgý hesaplandý ve distance degiþkenine atandý
                    distance = sqrt((S(i).xd - (S(n + 1).xd)) ^ 2 + (S(i).yd - (S(n + 1).yd)) ^ 2); %Öklit ile mesafe hesabý
                    C(cluster).distance = distance; %Cluster ýn distance kolonuna uzaklýgý atadýk.
                    C(cluster).id = i; %Clusterýn id sine node numarasýný atadýk
                    X(cluster) = S(i).xd; %Clusterýn x kordinatýný X dizisine ata
                    Y(cluster) = S(i).yd; %Clusterýn y kordinatýný Y dizisine ata
                    cluster = cluster + 1; %Clusterý 1 artýrdý
                 
                    %Sensör düðümleri k-bit verilerini ilettiðinde, enerji daðýtýmý þöyledir:
                 
                    % ClusterHeadin baz istasyonuna olan mesafesi ve do deðeri
                    % dikkate alýnarak ne kadar enerji harcamasý gerektigine karar
                    % veriliyor ve  ch in Enerjisni if satýrýna göre ve formüldeki þekilde güncelledik
                    distance;
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * (k) + Emp * k * (distance * distance * distance * distance));
                        S(i).E = S(i).E - ((ETX + EDA) * (k) + Emp * k * (distance * distance * distance * distance));
                    end
                    if (distance <= do)
                        S(i).E = S(i).E - ((ETX + EDA) * (k) + Efs * k * (distance * distance));
                        S(i).E = S(i).E - ((ETX + EDA) * (k) + Efs * k * (distance * distance));
                    end
                end
             
            end
        end
    end
 
    % Burada STATISTICS ve CLUSTERHS dizisine þuanki rounda oluþturulan clusterHead sayýsýný yazdýk
    % (1 eksiltmemizin sebebi cluster deðiþkeninin 1 den baþlamasý, dizi indeksi 1 den baþladýðý için)
    STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;
    CLUSTERHS(r + 1) = cluster - 1;
 
    %Normal nodlarýn hangi cluster head e baglanacagýný seçiyoruz.
 
    for i = 1:1:n
        if (S(i).type == 'N' && S(i).E > 0)
         
            %Aþagýdaki satýrlarda for döngüsünde baz istanyonu da dahil olmak üzere tüm clusterheadlerin elimizdeki noda olan uzaklýgýna bakarak mesafe ve node numarasýný tespit ettik
            %ve min_dis degiþkeninde mesafayi, min_dis_cluster deðiþkeninde ise bu clusterheadin numarasýný tuttuk
            if (cluster - 1 >= 1) %Cluster sayýsý birden fazlamý. cluster-1 ise gerçek cluster sayýsý normalde 0 tanýmlanamadýgý için
                min_dis = sqrt((S(i).xd - S(n + 1).xd) ^ 2 + (S(i).yd - S(n + 1).yd) ^ 2); %Elimizdeki node un baz istasyonuna olan uzaklýgý hesapalanarak mindis degiþkenine atandý
                min_dis_cluster = 1; %Degiþkenken oluþturmuþ 1 den baþlatmýþ
                for c = 1:1:cluster - 1 %Clusterhead kadar dön
                    temp = min(min_dis, sqrt((S(i).xd - C(c).xd) ^ 2 + (S(i).yd - C(c).yd) ^ 2)); %Burda node kendisine baz istasyonu mu yakýn yoksa var olan clusterheadlarden biri mi yakýn onu sorguluyo
                    if (temp < min_dis)
                        min_dis = temp;
                        min_dis_cluster = c; %Bize en yakýn olan clusterhead in numarasý minimum cluster degiþkeninde tut.
                    end
                end
             
                %Elimizdeki node az önce tespit ettiðimiz kendisine en yakýn Clusterheade paket göndererek enerji harcýyor
                min_dis; %min_dis degiþkeni do degerinden büyükse 1.if deki kadar enerji harcýyor
                if (min_dis > do)
                    S(i).E = S(i).E - (ETX * (k) + Emp * k * (min_dis * min_dis * min_dis * min_dis));
                end
                if (min_dis <= do); %min_dis degiþkeni do degerinden küyükse 2.if deki kadar enerji harcýyor
                    S(i).E = S(i).E - (ETX * (k) + Efs * k * (min_dis * min_dis));
                end
             
                %Cluster headlerin baz veri göndererek enerjilerinin güncellenmesi
                if (min_dis > 0)
                    %Elimizdeki clusterhead ile baz istasyonu arasýndaki mesafeyi hesapladýk
                    distance = sqrt((S(C(min_dis_cluster).id).xd - (S(n + 1).xd)) ^ 2 + (S(C(min_dis_cluster).id).yd - (S(n + 1).yd)) ^ 2);
                 
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * k); %ClusterHeadin yukarda ki node dan veri toplarken kaybettiði enerjiyi güncelledik
                 
                    %Clusterhead baz istasyonuna olan uzaklýgýna göre baz istanyonuna veri gönderip enerji harcýyor
                    if (distance > do)
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (k) + Emp * k * (distance * distance * distance * distance));
                    end
                    if (distance <= do)
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (k) + Efs * k * (distance * distance));
                    end
                 
                    PACKETS_TO_CH(r + 1) = n - total_dead - cluster + 1; %?
                end
             
                S(i).min_dis = min_dis; % node dizisinin,içeresinde elimizdeki nodun min_dis özelliðiðe kendisine en yakýn olan clusterýnhead in mesafesini atatýk
                S(i).min_dis_cluster = min_dis_cluster; % node dizisinin,içeresinde elimizdeki nodun min_dis_cluster özelliðiðe kendisine en yakýn olan clusterýnhead in id degerini atatýk
             
            end
        end
    end
    hold on;
 
    countCHs;
    rcountCHs = rcountCHs + countCHs;
 
    %Her bir roun için, elimizdeki nodelarýn enerjileri toplayýp node sayýsýna bölerek o roundun ortalama enerciyi hesapladýk
    sum = 0;
    for i = 1:1:n
        if (S(i).E > 0)
            sum = sum + S(i).E;
        end
    end
    avg = sum / n;
    STATISTICS(r + 1).AVG = avg; %STATISTICS degiþkeninde rounda ait ortalama enerjiyi yazdýk
    sum;
 
    %
    % Code for Voronoi Cells
    % Unfortynately if there is a small
    % number of cells, Matlab's voronoi
    % procedure has some problems
    %
    % [vx,vy]=voronoi(X,Y);
    % plot(X,Y,'r*',vx,vy,'b-');
    % hold on;
    %
    % voronoi(X,Y);
    % axis([0 xm 0 ym]);
 
end
plot(S(n + 1).xd, S(n + 1).yd, 'x'); %Baz istasyonunu ekrana x olarak çizdirdik

figure(3);
for r = 0:1:rmax - 1
    ylabel('Her Nodeun Ortalama Enerjisi');
    xlabel('Round Numarasý');
    plot([r r + 1], [STATISTICS(r + 1).AVG STATISTICS(r + 2).AVG], 'red');
    hold on;
end

figure(4);
for r = 0:1:rmax - 1
    ylabel('Ölü Node Sayýsý');
    xlabel('Round Numarasý');
    plot([r r + 1], [STATISTICS(r + 1).DEAD STATISTICS(r + 2).DEAD], 'red');
    hold on;
end